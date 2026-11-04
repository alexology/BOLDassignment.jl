using BOLDassignment
using DataFrames
using FASTX
using HTTP
using JSON3
using Test

include("mock.jl")
include("atomic_publication.jl")

const S1 = "ACGTACGTAAACCCGGGTTT"
const S2 = "TTTTGGGGCCCCAAAATTGA"
const S3 = "GGGAAATTTCCCGGGAAATT"
const S4 = "CATCATCATCATCATCATCA"
const S5 = "AAAAAAAAAACCCCCCCCCC"

catalog() = Dict{String,Vector{Any}}(
    S1 => [("P1", "BOLD:A1", 0.0, beetle("Onthophagus", "Onthophagus taurus")),
           ("P2", "BOLD:A1", 1.0, beetle("Onthophagus", "Onthophagus taurus")),
           ("P3", "BOLD:A2", 2.5, beetle("Onthophagus", "Onthophagus illyricus"))],
    S2 => [("P4", "BOLD:B1", 6.0, beetle("Aphodius", "Aphodius fimetarius"))],
    S3 => [],
    S4 => [("P5", "BOLD:C1", 0.5, beetle("Copris", "Copris lunaris"))],
    S5 => [("P6", "BOLD:D1", 8.0, lineage("Arthropoda", "Insecta", "Diptera", "Muscidae", "", ""))],
)

@testset "BOLDassignment" begin
    @testset "Sequence normalization" begin
        @test normalize_sequence(" ac-g.u\n") == "ACGT"
        @test normalize_sequence("nRy") == "NRY"
        @test_throws ArgumentError normalize_sequence("ACGZ")
        @test_throws ArgumentError normalize_sequence("--")
    end

    @testset "Valid names" begin
        for good in ["Onthophagus taurus", "Onthophagus (Palaeonthophagus) taurus",
                     "Nymphalis c-album", "Carabus violaceus purpurascens"]
            @test is_valid_name(:species, good)
        end
        for bad in ["Onthophagus sp.", "Onthophagus sp", "Onthophagus spp.", "Onthophagus sp. 1",
                    "Onthophagus cf. taurus", "Onthophagus cf taurus", "Onthophagus aff. taurus",
                    "Onthophagus nr. taurus", "Onthophagus taurus complex", "Diptera sp. BIOUG1234-A01",
                    "Cecidomyiidae sp. 4 LH-2012", "Onthophagus", "Onthophagus Taurus",
                    "Onthophagus taurus1", "Onthophagus 'taurus'", "", missing]
            @test !is_valid_name(:species, bad)
        end
        @test is_valid_name(:genus, "Onthophagus")
        @test is_valid_name(:family, "Scarabaeidae")
        for bad in ["Unknown", "Malaise4919", "Onthophagus sp.", "incertae sedis", "", missing]
            @test !is_valid_name(:genus, bad)
        end
    end

    @testset "LCA with agreement" begin
        lin(g, s) = (phylum = "Arthropoda", class = "Insecta", order = "Coleoptera",
                      family = "Scarabaeidae", genus = g, species = s)
        sims = [100.0, 99.0, 98.0]
        lins = [lin("Alpha", "Alpha alba"), lin("Alpha", "Alpha alba"), lin("Alpha", "Alpha beta")]
        path = lca_path(sims, lins)
        @test last(path).rank == :genus && last(path).taxon == "Alpha"
        @test all(p -> p.agreement == 100, path)
        path = lca_path(sims, lins; min_agreement = 66)
        @test last(path).rank == :species && last(path).taxon == "Alpha alba"
        @test last(path).agreement ≈ 200 / 3
        @test last(path).n_support == 2 && last(path).n_pool == 3
        # 2/3 reaches exactly min_agreement = 200/3.
        @test last(lca_path(sims, lins; min_agreement = 200 / 3)).rank == :species
        # The 98% hit is out of the pool with top_delta = 1.
        @test last(lca_path(sims, lins; top_delta = 1.0)).taxon == "Alpha alba"
        # Species threshold 97: a 96% hit does not take part at species rank.
        alba_beta = [lin("Alpha", "Alpha alba"), lin("Alpha", "Alpha beta")]
        @test last(lca_path([100.0, 96.0], alba_beta)).taxon == "Alpha alba"
        # Ties stop the assignment even when they reach min_agreement.
        path = lca_path([100.0, 100.0], alba_beta; min_agreement = 50)
        @test last(path).rank == :genus
        # Hits not identified at a rank leave the pool, unless ignore_missing = false.
        lins = [lin("Alpha", "Alpha alba"), lin("Alpha", missing)]
        @test last(lca_path([100.0, 100.0], lins)).taxon == "Alpha alba"
        @test last(lca_path([100.0, 100.0], lins; ignore_missing = false)).rank == :genus
        # Open nomenclature above the species threshold: no vote at species, vote at genus.
        lins = [lin("Onthophagus", "Onthophagus taurus"), lin("Onthophagus", "Onthophagus sp."),
                lin("Onthophagus", "Onthophagus cf. illyricus"),
                lin("Onthophagus", "Onthophagus sp. BIOUG123-A01")]
        path = lca_path([99.0, 99.5, 99.2, 98.0], lins)
        @test last(path).taxon == "Onthophagus taurus"
        @test last(path).n_pool == 1
        @test path[end-1].rank == :genus && path[end-1].n_pool == 4
        @test last(lca_path([99.0, 99.5, 99.2, 98.0], lins; ignore_missing = false)).rank == :genus
        # Only open-nomenclature hits: assignment stops at genus.
        @test last(lca_path([99.5, 99.2], lins[2:3])).rank == :genus
        # A hit identified only to family does not vote at genus nor at species.
        lins = [lin("Alpha", "Alpha alba"), lin(missing, "Alpha beta")]
        path = lca_path([99.0, 99.0], lins)
        @test last(path).taxon == "Alpha alba" && last(path).n_pool == 1
        # Custom name check.
        @test last(lca_path([99.0, 99.0], [lin("Alpha", "Alpha alba"), lin("Alpha", "Alpha sp.")];
                            valid_name = (r, n) -> !ismissing(n))).rank == :genus
        # Below the phylum threshold nothing is assigned.
        @test isempty(lca_path([40.0], [lin("Alpha", "Alpha alba")]))
        @test isempty(lca_path(Float64[], typeof(lin("Alpha", "Alpha alba"))[]))
        # A deeper rank must agree with the ranks already assigned.
        full(o, f, g, s) = (phylum = "Arthropoda", class = "Insecta", order = o, family = f,
                            genus = g, species = s)
        lins = [full("Coleoptera", "Carabidae", "Carabus", "Carabus alpha"),
                full("Diptera", "Muscidae", "Carabus", "Carabus alpha"),
                full("Coleoptera", "Carabidae", "Nebria", "Nebria beta")]
        path = lca_path([99.0, 99.0, 99.0], lins; min_agreement = 60)
        @test [p.rank for p in path] == [:phylum, :class, :order, :family]
        @test_throws ArgumentError lca_path([1.0], lins[1:1]; min_agreement = 0)
        @test_throws ArgumentError lca_path([1.0], lins; min_agreement = 50)
        # Thresholds given for some ranks keep the defaults for the others.
        alba = [lin("Alpha", "Alpha alba")]
        @test last(lca_path([96.0], alba; thresholds = Dict(:species => 96.0))).rank == :species
        @test last(lca_path([96.0], alba; thresholds = Dict(:species => 99.0))).rank == :genus
        @test last(lca_path([94.0], alba; thresholds = Dict(:species => 90.0))).rank == :family
        @test_throws ArgumentError lca_path([99.0], alba; thresholds = Dict(:specie => 98.0))
        @test_throws ArgumentError lca_path([99.0], alba; thresholds = Dict(:species => NaN))
        @test_throws ArgumentError lca_path([99.0], alba; thresholds = Dict(:species => 101))
        @test_throws ArgumentError lca_path([NaN], alba)
        @test_throws ArgumentError lca_path([100.5], alba)
        @test_throws ArgumentError lca_path([99.0], alba; top_delta = NaN)
        @test_throws ArgumentError lca_path([99.0], alba; ranks = (:genus, :genus))
        @test_throws ArgumentError lca_path([99.0], alba; ranks = (:genus, :taxon))
        @test_throws ArgumentError lca_path([99.0], alba; ranks = (:genus, :subgenus))
        sub = [(genus = "Alpha", subgenus = "Beta", species = "Alpha alba")]
        path = lca_path([99.0], sub; ranks = (:genus, :subgenus, :species),
                        thresholds = Dict(:subgenus => 96.0))
        @test [p.rank for p in path] == [:genus, :subgenus, :species]
    end

    @testset "Identification with a fake engine" begin
        fake = FakeEngine(catalog(); polls_before_done = 2)
        engine = fake_engine(fake)
        hits = bold_identify([S1, lowercase(S2), S3, S1, missing, S4, S5]; chunk_size = 2, engine)
        @test length(fake.submissions) == 3
        @test sort(unique(hits.sequence)) == sort([S1, S2, S3, S4, S5])
        @test nrow(hits) == 3 + 1 + 1 + 1 + 1
        s1 = hits[hits.sequence .== S1, :]
        @test sort(s1.similarity) == [97.5, 99.0, 100.0]
        @test Set(s1.bin_uri) == Set(["BOLD:A1", "BOLD:A2"])
        @test all(==("public"), s1.record_status)
        none = hits[hits.sequence .== S3, :]
        @test ismissing(only(none.process_id)) && ismissing(only(none.similarity))
        @test ismissing(only(hits[hits.sequence .== S5, :species]))
        # Request spacing respects the documented limits.
        posts = [t for (t, m, _) in fake.log if m == "POST"]
        @test all(diff(posts) .>= engine.submit_interval)
        times = [t for (t, _, _) in fake.log]
        @test all(diff(times) .>= engine.request_interval)
        @test_throws ArgumentError bold_identify([S1]; chunk_size = 1001, engine)
        @test_throws ArgumentError bold_identify([S1]; db = :nope, engine)
        @test nrow(bold_identify(String[]; engine)) == 0
    end

    @testset "Retries and checkpoints" begin
        fake = FakeEngine(catalog())
        fake.fail_next = [503, 429]
        dir = mktempdir()
        hits = bold_identify([S1, S2, S4]; chunk_size = 2, checkpoint = dir, engine = fake_engine(fake))
        @test length(readdir(dir)) == 3
        @test all(endswith(".jsonl"), readdir(dir))
        n_requests = length(fake.log)
        # The checkpoint does not depend on the input order: everything is reused.
        same = bold_identify([S4, lowercase(S1), S2, S4]; chunk_size = 2, checkpoint = dir,
                             engine = fake_engine(fake))
        @test length(fake.log) == n_requests
        @test isequal(same, hits)
        # Adding a sequence submits only that sequence; the others are reused.
        n_subs = length(fake.submissions)
        more = bold_identify([S1, S3, S2, S4]; chunk_size = 2, checkpoint = dir, engine = fake_engine(fake))
        @test length(fake.submissions) == n_subs + 1
        @test fake.submissions["sub$(n_subs + 1)"] == ["1" => S3]
        @test isequal(more[more.sequence .!= S3, :], hits)
        @test length(readdir(dir)) == 4
        key = BOLDassignment._sequence_key(BOLDassignment._database_code(:animal_public),
                                           BOLDassignment._mode(:rapid), S3)
        rm(joinpath(dir, key * ".jsonl"))
        # A corrupted checkpoint is reported, resubmitted and rewritten.
        files = sort(readdir(dir; join = true))
        write(files[1], "{\"seqid\":\"1\",\"res")
        fixed = @test_logs (:warn, r"Unreadable checkpoint") match_mode = :any bold_identify(
            [S1, S2, S4]; chunk_size = 2, checkpoint = dir, engine = fake_engine(fake))
        @test isequal(fixed, hits)
        @test length(fake.log) > n_requests
        @test sort(readdir(dir; join = true)) == files
        @test !any(endswith(".partial"), readdir(dir))
        @test_throws ArgumentError bold_identify([S1]; checkpoint = files[1],
                                                 engine = fake_engine(fake))
        # A client error other than 408/429 is not retried.
        bad = FakeEngine(catalog())
        bad.fail_next = [400]
        @test_throws ErrorException bold_identify([S1]; engine = fake_engine(bad))
        @test length(bad.log) == 1
        # Submissions that never complete are resubmitted a bounded number of times.
        stuck = FakeEngine(catalog(); polls_before_done = typemax(Int))
        @test_throws ErrorException bold_identify([S1]; engine = fake_engine(stuck; timeout = 100.0,
                                                                            max_resubmits = 2))
        @test length(stuck.submissions) == 3
        err = try
            bold_identify([S1]; checkpoint = mktempdir(),
                          engine = fake_engine(stuck; timeout = 100.0, max_resubmits = 0))
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("last submission sub4", err.msg) && occursin("call again to resume", err.msg)
        # A finished submission without the line of a sequence: reported once,
        # returned without hits and not checkpointed, so it is submitted again.
        partial = FakeEngine(catalog())
        partial.omit = [S2]
        pdir = mktempdir()
        got = @test_logs (:warn, r"without the result of 1 sequence") match_mode = :any bold_identify(
            [S1, S2]; checkpoint = pdir, engine = fake_engine(partial))
        @test length(partial.submissions) == 1
        @test ismissing(only(got[got.sequence .== S2, :process_id]))
        @test count(==(S1), got.sequence) == 3
        @test only(got[got.sequence .== S2, :result_received]) == false
        @test all(got[got.sequence .== S1, :result_received])
        @test length(readdir(pdir)) == 1
        empty!(partial.omit)
        again = bold_identify([S1, S2]; checkpoint = pdir, engine = fake_engine(partial))
        @test partial.submissions["sub2"] == ["1" => S2]
        @test only(again[again.sequence .== S2, :process_id]) == "P4"
        # The spacing is shared by calls using the same throttle.
        fake = FakeEngine(catalog())
        shared = BOLDassignment._Throttle()
        bold_identify([S1]; engine = fake_engine(fake; throttle = shared))
        bold_identify([S2]; engine = fake_engine(fake; throttle = shared))
        posts = [t for (t, m, _) in fake.log if m == "POST"]
        @test length(posts) == 2 && posts[2] - posts[1] >= 21.0
        @test_throws ArgumentError IDEngine(request_interval = -1.0)
        @test_throws ArgumentError IDEngine(max_active = 0)
        @test IDEngine(base_url = "https://x.test/").base_url == "https://x.test"
    end

    @testset "Parsing of BOLD results" begin
        parse_hits(seqs, lines...) =
            BOLDassignment._hits_table(seqs, Vector{UInt8}(codeunits(join(lines, "\n"))))
        hit(pdist) = "{\"P1|COI-5P|BOLD:A1|x|public\":{\"pdist\":$pdist," *
                     "\"taxonomy\":{\"genus\":\"Copris\"}}}"
        line(k, results) = "{\"seqid\":\"$k\",\"results\":$results}"
        @test parse_hits([S1], line(1, hit(1.5))).similarity == [98.5]
        # An explicit empty result is received (no match); an absent line is not.
        @test parse_hits([S1], line(1, "{}")).result_received == [true]
        @test BOLDassignment._hits_table([S1, S2], [line(1, "{}"), nothing]).result_received == [true, false]
        @test_throws ErrorException parse_hits([S1, S2], line(1, "{}"), line(1, "{}"))
        @test_throws ErrorException parse_hits([S1, S2], line(1, "{}"))
        @test_throws ErrorException parse_hits([S1], line(3, "{}"))
        for bad in ("null", "-1", "101", "\"x\"")
            @test_throws ErrorException parse_hits([S1], line(1, hit(bad)))
        end
    end

    @testset "Assignment of a DataFrame" begin
        fake = FakeEngine(catalog())
        asvs = DataFrame(id = ["a", "b", "c", "d", "e"], ASV = [S1, S2, missing, S1, S3])
        tax = bold_assign(asvs; engine = fake_engine(fake))
        @test tax.id == asvs.id && isequal(tax.ASV, asvs.ASV)
        @test isequal(tax.assigned_rank, ["genus", "family", missing, "genus", missing])
        @test isequal(tax.taxon, ["Onthophagus", "Scarabaeidae", missing, "Onthophagus", missing])
        @test isequal(tax.n_hits, [3, 1, missing, 3, 0])
        @test tax.best_process_id[1] == "P1"
        @test tax.genus_agreement[1] == 100
        @test ismissing(tax.species_agreement[1])
        # Majority LCA reuses the hits without network access.
        hits = bold_identify(asvs.ASV; engine = fake_engine(fake))
        n_requests = length(fake.log)
        tax = bold_assign(asvs; hits, min_agreement = 60)
        majority = tax
        @test length(fake.log) == n_requests
        @test tax.species[1] == "Onthophagus taurus"
        @test tax.species_agreement[1] ≈ 200 / 3
        @test_throws ArgumentError bold_assign(rename(asvs, :ASV => :seq); hits)
        tax = bold_assign(rename(asvs, :ASV => :seq); sequence_col = :seq, hits, prefix = "bold_")
        @test hasproperty(tax, :bold_species)
        @test_throws ArgumentError bold_assign(hcat(asvs, DataFrame(taxon = 1:5)); hits)
        @test_throws ArgumentError bold_assign(asvs; hits, min_agrement = 60)
        @test_throws ArgumentError bold_assign(DataFrame(ASV = ["GGGG"]); hits)
        # Repeated hit rows are counted once; conflicting repeats are refused.
        tax2 = @test_logs (:warn, r"Ignored 4 repeated") bold_assign(asvs; hits = vcat(hits, hits),
                                                                     min_agreement = 60)
        @test isequal(tax2, majority)
        conflict = copy(hits)
        conflict[findfirst(isequal("P1"), conflict.process_id), :similarity] = 90.0
        @test_throws ArgumentError bold_assign(asvs; hits = vcat(hits, conflict))
        # Hits with sequences written differently still match.
        lower = transform(hits, :sequence => ByRow(lowercase) => :sequence)
        @test isequal(bold_assign(asvs; hits = lower, min_agreement = 60), majority)
        # Arguments are checked before any request.
        silent = FakeEngine(catalog())
        engine = fake_engine(silent)
        @test_throws ArgumentError bold_assign(asvs; engine, thresholds = Dict(:specie => 98))
        @test_throws ArgumentError bold_assign(asvs; engine, min_agreement = 0)
        @test_throws ArgumentError bold_assign(hcat(asvs, DataFrame(species_agreement = 1:5)); engine)
        @test_throws r"row 2 of :ASV" bold_assign(DataFrame(ASV = [S1, "ACGZ"]); engine)
        @test isempty(silent.log)
        empty_tax = bold_assign(DataFrame(ASV = [missing, missing]); engine)
        @test isempty(silent.log)
        @test all(ismissing, empty_tax.taxon) && hasproperty(empty_tax, :species_agreement)
    end

    @testset "Real BOLD response" begin
        # Response of the identification engine for CASBI15164-26 (Onthophagus
        # taurus) and a low-complexity sequence, public.tax-derep, rapid mode,
        # retrieved on 2026-10-23; seqids renumbered as the package submits them.
        body = read(joinpath(@__DIR__, "data", "bold_results_taurus.jsonl"))
        seqs = [normalize_sequence(String(JSON3.read(l)[:sequence])) for l in eachline(IOBuffer(body))]
        hits = BOLDassignment._hits_table(seqs, body)
        @test nrow(hits) == 15
        real = hits[.!ismissing.(hits.similarity), :]
        @test nrow(real) == 14 && all(==("BOLD:AAH6811"), real.bin_uri)
        @test all(==("public"), real.record_status)
        @test count(==(100.0), real.similarity) == 9
        @test real.similarity[real.process_id .== "TDAOE1457-22"] ≈ [100 - 0.3184713375796178]
        @test count(ismissing, real.species) == 4
        @test all(isequal("Onthophagus"), real.genus)
        # Alignment fields kept from the raw hits.
        @test all(==("COI-5P"), real.marker) && !any(ismissing, real.taxid)
        @test real.taxid[real.process_id .== "UTOCP034-25"] == ["153867"]
        @test real.query_coverage[real.process_id .== "GBCL3727-08"] == [91.0]
        @test real.align_length[real.process_id .== "UTOCP034-25"] == [314]
        @test real.query_coverage[real.process_id .== "UTOCP034-25"] == [100.0]
        @test all(x -> 0 < x <= 100, real.query_coverage)
        @test all(ismissing, hits.marker[ismissing.(hits.similarity)])
        @test eltype(hits.align_length) == Union{Missing,Int}
        # Strict LCA: one O. illyricus record at 99.68 % keeps the species open.
        strict = lca_assign(hits)
        @test isequal(strict.assigned_rank, ["genus", missing])
        @test strict.genus_agreement[1] == 100
        # 9 of the 10 records identified to species are O. taurus; the four
        # records without species do not vote (BOLD classification: 90.0).
        majority = lca_assign(hits; min_agreement = 90)
        @test majority.species[1] == "Onthophagus taurus"
        @test majority.species_agreement[1] == 90
        @test only(lca_path(Float64.(real.similarity),
                            [NamedTuple{RANKS}(Tuple(r[k] for k in RANKS)) for r in eachrow(real)];
                            min_agreement = 90)[end:end]).n_pool == 10
        @test lca_assign(hits; min_agreement = 90, ignore_missing = false).assigned_rank[1] == "genus"
        @test majority.n_hits == [14, 0]
    end

    @testset "Status formats" begin
        reply(body) = IDEngine(; request = (m, u, p) -> (200, Vector{UInt8}(codeunits(body))),
                               request_interval = 0.0, verbose = false,
                               throttle = BOLDassignment._Throttle())
        finished(body) = BOLDassignment._finished(reply(body), "id")
        @test finished("{\"status\":\"complete\"}")
        @test !finished("{\"status\":\"queued\"}")
        @test !finished("{\"queued\":0,\"processing\":2,\"completed\":0}")
        @test finished("{\"queued\":0,\"processing\":0,\"completed\":2}")
        @test !finished("{}")
        @test !finished("not json")
    end

    @testset "Portal batches" begin
        ids = ["PROC$(lpad(i, 5, '0'))-26" for i in 1:40]
        batches = BOLDassignment._portal_batches(ids)
        @test reduce(vcat, batches) == ids
        @test all(b -> length(join(("ids:processid:" * p for p in b), ';')) <= 250, batches)
        @test length(batches) > 1
        # An ID that fits alone is accepted after another one (236 + 14 = 250).
        long = "L" ^ 236
        @test BOLDassignment._portal_batches(["P1", long]) == [["P1"], [long]]
        @test_throws ArgumentError BOLDassignment._portal_batches(["L" ^ 237])
    end

    @testset "Portal records with a fake portal" begin
        # Portal shaped as observed on 2026-10-23: /api/query returns a query_id,
        # /api/documents/{id} pages of {"data": [...]}.
        db = Dict("P$i" => Dict("processid" => "P$i", "species" => i == 2 ? nothing : "Sp $i",
                                "coord" => [41.988, -120.036]) for i in 1:250)
        urls = String[]
        function portal(method, url, payload)
            push!(urls, url)
            if occursin("/api/query?", url)
                q = HTTP.unescapeuri(match(r"query=([^&]*)", url)[1])
                ids = [last(split(t, ':')) for t in split(q, ';')]
                return (200, Vector{UInt8}(codeunits(JSON3.write(Dict("query_id" => join(ids, ','))))))
            end
            m = match(r"/api/documents/([^?]*)\?length=(\d+)&start=(\d+)", url)
            ids = split(HTTP.unescapeuri(m[1]), ',')
            found = [db[i] for i in ids if haskey(db, i)]
            n, start = parse(Int, m[2]), parse(Int, m[3])
            page = found[start+1:min(end, start + n)]
            return (200, Vector{UInt8}(codeunits(JSON3.write(Dict("data" => page)))))
        end
        engine = IDEngine(; request = portal, request_interval = 0.0, verbose = false,
                          throttle = BOLDassignment._Throttle())
        recs = bold_records(["P1", "P2", "missingID", missing, " P1 "]; engine,
                            base_url = "https://portal.test")
        @test recs.processid == ["P1", "P2"]
        @test isequal(recs.species, ["Sp 1", missing])
        @test eltype(recs.processid) == String
        @test eltype(recs.species) == Union{Missing,String}
        @test recs.coord == fill("[41.988,-120.036]", 2)
        @test all(startswith("https://portal.test/api/"), urls)
        # Many process IDs are split in several queries (250-character limit).
        many = bold_records(["P$i" for i in 1:250]; engine, base_url = "https://portal.test",
                            fields = [:processid])
        @test many.processid == ["P$i" for i in 1:250]
        @test propertynames(many) == [:processid]
        @test isempty(bold_records(String[]; engine))
        # One document per marker: the `marker` keyword keeps only the chosen one.
        db["P1"]["marker_code"] = "COI-5P"
        db["P2"]["marker_code"] = "COI-3P"
        coi5 = bold_records(["P1", "P2", "P3"]; engine, base_url = "https://portal.test",
                            marker = "coi-5p")
        @test coi5.processid == ["P1"]
        @test nrow(bold_records(["P1", "P2"]; engine, base_url = "https://portal.test")) == 2
    end

    @testset "Query coverage fallbacks" begin
        rec(s) = JSON3.read(s)
        @test BOLDassignment._query_coverage(rec("{\"qcovhsp\": 87}"), 300) == 87.0
        @test BOLDassignment._query_coverage(rec("{\"qstart\": 1, \"qend\": 150}"), 300) == 50.0
        @test BOLDassignment._query_coverage(rec("{\"qstart\": 300, \"qend\": 1}"), 300) == 100.0
        @test ismissing(BOLDassignment._query_coverage(rec("{\"qcovhsp\": 120}"), 300))
        @test ismissing(BOLDassignment._query_coverage(rec("{}"), 300))
        @test BOLDassignment._align_length(rec("{\"length\": 314}")) == 314
        @test ismissing(BOLDassignment._align_length(rec("{\"length\": 0}")))
        @test ismissing(BOLDassignment._align_length(rec("{\"length\": \"x\"}")))
    end

    if get(ENV, "BOLDASSIGNMENT_LIVE", "") == "1"
        @testset "Live BOLD request" begin
            include("live.jl")
        end
    end
end
