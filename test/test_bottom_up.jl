# Bottom-up direction of the LCA and invariance of the default direction.

include("legacy_lca.jl")

# Deterministic pseudo-random numbers, independent of the Julia version.
mutable struct TestRNG
    state::UInt64
end
function next!(g::TestRNG)
    g.state = g.state * 6364136223846793005 + 1442695040888963407
    return g.state >> 33
end
uniform!(g::TestRNG) = next!(g) / 2.0^31                  # in [0, 1)
pick!(g::TestRNG, xs) = xs[Int(next!(g) % length(xs)) + 1]
chance!(g::TestRNG, p) = uniform!(g) < p

const BU_RANKS = BOLDassignment.RANKS

lineage6(p, c, o, f, g, s) = (phylum = p, class = c, order = o, family = f, genus = g, species = s)

# Random hit sets over a small taxonomy with open nomenclature, missing names
# and inconsistent classifications (the same species under two genera).
function random_hits!(rng::TestRNG)
    n = Int(next!(rng) % 31)
    sims = Float64[]
    lins = []
    for _ in 1:n
        p, c, o, f = pick!(rng, ("Arthropoda", "Mollusca")), pick!(rng, ("Insecta", "Clitellata")),
                     pick!(rng, ("Plecoptera", "Diptera")), pick!(rng, ("Nemouridae", "Leuctridae"))
        g = pick!(rng, ("Protonemura", "Nemoura", "Leuctra"))
        epithet = pick!(rng, ("intricata", "cinerea", "fusca"))
        s = string(g, " ", epithet)
        chance!(rng, 0.08) && (g = pick!(rng, ("Protonemura", "Nemoura", "Leuctra")))  # inconsistent
        chance!(rng, 0.10) && (s = pick!(rng, (string(g, " sp."), string(g, " cf. ", epithet),
                                                string(g, " aff. ", epithet),
                                                string(g, " sp. BIOUG1234-A01"), missing)))
        chance!(rng, 0.05) && (g = pick!(rng, (missing, "Unknown")))
        chance!(rng, 0.03) && (f = missing)
        chance!(rng, 0.02) && (o = missing)
        sim = chance!(rng, 0.2) ? pick!(rng, (97.0, 95.0, 90.0, 85.0, 100.0)) :
              round(80 + 20 * uniform!(rng); digits = 1)
        push!(sims, sim)
        push!(lins, lineage6(p, c, o, f, g, s))
    end
    return sims, lins
end

function random_options!(rng::TestRNG)
    thresholds = chance!(rng, 0.7) ? DEFAULT_THRESHOLDS :
        Dict(:species => pick!(rng, (96.0, 97.0, 98.5)), :genus => pick!(rng, (93.0, 95.0)))
    return (min_agreement = pick!(rng, (100, 95, 90, 80, 200 / 3, 60, 50.5)),
            top_delta = pick!(rng, (Inf, Inf, 0.0, 1.0, 2.5, 5.0)),
            ignore_missing = chance!(rng, 0.7), thresholds = thresholds)
end

# Independent reference implementation of the bottom-up rule, written
# directly from its description.
function oracle_bottom_up(sims, lins, opts)
    limits = BOLDassignment._check_thresholds(opts.thresholds, BU_RANKS)
    valid = BOLDassignment.is_valid_name
    out = NamedTuple{(:rank, :taxon, :agreement, :n_support, :n_pool),
                     Tuple{Symbol,String,Float64,Int,Int}}[]
    isempty(sims) && return out, false
    floor_sim = maximum(sims) - opts.top_delta
    for j in length(BU_RANKS):-1:1
        r = BU_RANKS[j]
        inpool = [s >= limits[r] && s >= floor_sim for s in sims]
        named = [valid(r, l[r]) for l in lins]
        votes = [lins[i][r] for i in eachindex(sims) if inpool[i] && named[i]]
        isempty(votes) && continue
        npool = opts.ignore_missing ? length(votes) : count(inpool)
        names = unique(votes)
        counts = [count(==(x), votes) for x in names]
        m = maximum(counts)
        count(==(m), counts) == 1 || continue
        100 * m / npool >= opts.min_agreement * (1 - 1e-12) || continue
        taxon = names[argmax(counts)]
        support = [i for i in eachindex(sims) if inpool[i] && named[i] && lins[i][r] == taxon]
        push!(out, (rank = r, taxon = taxon, agreement = 100 * m / npool, n_support = m,
                    n_pool = npool))
        conflict = false
        for q in j-1:-1:1
            rq = BU_RANKS[q]
            v = [lins[i][rq] for i in support if valid(rq, lins[i][rq])]
            if isempty(v)
                conflict = true
                break
            end
            nm = unique(v)
            cc = [count(==(x), v) for x in nm]
            mm = maximum(cc)
            if count(==(mm), cc) > 1
                conflict = true
                break
            end
            mm < length(support) && (conflict = true)
            pushfirst!(out, (rank = rq, taxon = nm[argmax(cc)],
                             agreement = 100 * mm / length(support), n_support = mm,
                             n_pool = length(support)))
        end
        return out, conflict
    end
    return out, false
end

# Distinct DNA strings for the sequence column of hit tables.
seq_of(k) = join("ACGT"[((k >> (2i)) & 3) + 1] for i in 0:11)

function hit_table(sets)
    rows = DataFrame(sequence = String[], similarity = Union{Missing,Float64}[],
                     process_id = Union{Missing,String}[])
    for r in BU_RANKS
        rows[!, r] = Union{Missing,String}[]
    end
    for (k, (sims, lins)) in enumerate(sets)
        if isempty(sims)
            push!(rows, merge((sequence = seq_of(k), similarity = missing, process_id = missing),
                              NamedTuple{BU_RANKS}(ntuple(_ -> missing, length(BU_RANKS)))))
        end
        for (i, (s, l)) in enumerate(zip(sims, lins))
            push!(rows, merge((sequence = seq_of(k), similarity = s, process_id = "P$(k)_$i"), l))
        end
    end
    return rows
end

ins(g, s) = lineage6("Arthropoda", "Insecta", "Plecoptera", "Nemouridae", g, s)
bu(sims, lins; kw...) = lca_path(sims, lins; direction = :bottom_up, kw...)

@testset "Bottom-up LCA: examples" begin
    # Protonemura: 10 P. intricata 97.0–99.3 %, 35 other Plecoptera 85.9–96.9 %,
    # 55 Diptera 85.4–86.5 %.
    sims = Float64[]
    lins = []
    for k in 0:9
        push!(sims, 97.0 + 2.3 * (9 - k) / 9); push!(lins, ins("Protonemura", "Protonemura intricata"))
    end
    for k in 0:34
        push!(sims, 85.9 + 11.0 * k / 34)
        push!(lins, lineage6("Arthropoda", "Insecta", "Plecoptera", k % 2 == 0 ? "Nemouridae" : "Leuctridae",
                             k % 2 == 0 ? "Nemoura" : "Leuctra",
                             k % 2 == 0 ? "Nemoura cinerea" : "Leuctra fusca"))
    end
    for k in 0:54
        push!(sims, 85.4 + 1.1 * k / 54)
        push!(lins, lineage6("Arthropoda", "Insecta", "Diptera", "Chironomidae", "Chironomus",
                             "Chironomus riparius"))
    end
    @test length(sims) == 100 && maximum(sims) ≈ 99.3
    @test count(l -> l.order == "Plecoptera", lins) == 45
    @test all(collect(extrema(sims[11:45])) .≈ [85.9, 96.9]) && all(s -> 85.4 <= s <= 86.5 + 1e-9, sims[46:end])
    top = lca_path(sims, lins; top_delta = Inf)
    @test last(top).rank == :class && last(top).taxon == "Insecta"
    @test lca_path(sims, lins) == top                       # default direction
    path = bu(sims, lins)
    @test last(path).rank == :species && last(path).taxon == "Protonemura intricata"
    @test (last(path).n_support, last(path).n_pool, last(path).agreement) == (10, 10, 100.0)
    @test [p.taxon for p in path] == ["Arthropoda", "Insecta", "Plecoptera", "Nemouridae",
                                      "Protonemura", "Protonemura intricata"]
    @test all(p -> p.n_pool == 10 && p.agreement == 100, path)
    t = lca_assign(hit_table([(sims, lins)]); direction = :bottom_up)
    @test (t.assigned_rank[1], t.decision_rank[1], t.lineage_conflict[1]) == ("species", "species", false)
    @test t.order[1] == "Plecoptera" && t.order_agreement[1] == 100
    @test lca_assign(hit_table([(sims, lins)])).assigned_rank[1] == "class"

    # A valid congener at 97.5 % blocks the species in both directions.
    two = [ins("Protonemura", "Protonemura intricata"), ins("Protonemura", "Protonemura nimborum")]
    @test last(lca_path([100.0, 97.5], two)).rank == :genus
    path = bu([100.0, 97.5], two)
    @test last(path).rank == :genus && last(path).taxon == "Protonemura"
    @test (last(path).n_support, last(path).n_pool) == (2, 2)

    # min_agreement = 95: 19 + 1 gives the species; 9 + 1 does not.
    mk(n) = vcat(fill(ins("Protonemura", "Protonemura intricata"), n),
                 [ins("Protonemura", "Protonemura nimborum")])
    sims20 = [97.0 + 3 * k / 19 for k in 0:19]
    for direction in (:top_down, :bottom_up)
        p = lca_path(sims20, mk(19); min_agreement = 95, direction)
        @test last(p).taxon == "Protonemura intricata" && last(p).agreement == 95
        p = lca_path(sims20[1:10], mk(9); min_agreement = 95, direction)
        @test last(p).rank == :genus
    end
    @test last(bu(sims20[1:10], mk(9); min_agreement = 95)).n_pool == 10

    # Open nomenclature and interim codes above 97 % do not vote for the
    # species but vote for the genus.
    opennames = [ins("Alpha", "Alpha alba"), ins("Alpha", "Alpha alba"), ins("Alpha", "Alpha alba"),
            ins("Alpha", "Alpha sp."), ins("Alpha", "Alpha cf. alba"), ins("Alpha", "Alpha aff. alba"),
            ins("Alpha", "Alpha sp. BIOUG1234-A01")]
    s7 = fill(98.0, 7)
    p = bu(s7, opennames)
    @test last(p).taxon == "Alpha alba" && (last(p).n_support, last(p).n_pool) == (3, 3)
    p = bu(s7, opennames; ignore_missing = false)                 # 3/7 at the species
    @test last(p).rank == :genus && (last(p).n_support, last(p).n_pool) == (7, 7)
    tie = vcat(opennames, fill(ins("Alpha", "Alpha beta"), 3))
    p = bu(fill(98.0, 10), tie)                               # species tie 3–3
    @test last(p).rank == :genus && (last(p).n_support, last(p).n_pool) == (10, 10)
    @test last(lca_path(fill(98.0, 10), tie)).rank == :genus  # top-down agrees here

    # Tie at the species: up to the genus.
    p = bu(fill(99.0, 4), [ins("Alpha", "Alpha alba"), ins("Alpha", "Alpha alba"),
                           ins("Alpha", "Alpha beta"), ins("Alpha", "Alpha beta")])
    @test last(p).rank == :genus

    # Empty species pool with a valid genus.
    p = bu([96.5, 95.2], [ins("Alpha", "Alpha alba"), ins("Alpha", "Alpha beta")])
    @test last(p).rank == :genus && last(p).taxon == "Alpha"

    # Supporting hits under two genera: majority with conflict, or stop on a tie.
    mixed = [ins("Alpha", "Alpha alba"), ins("Alpha", "Alpha alba"), ins("Alpha", "Alpha alba"),
             ins("Beta", "Alpha alba")]
    p = bu(fill(99.0, 4), mixed)
    @test last(p).taxon == "Alpha alba"
    g = only(filter(x -> x.rank == :genus, p))
    @test (g.taxon, g.n_support, g.n_pool, g.agreement) == ("Alpha", 3, 4, 75.0)
    @test only(filter(x -> x.rank == :family, p)).agreement == 100
    t = lca_assign(hit_table([(fill(99.0, 4), mixed)]); direction = :bottom_up)
    @test t.lineage_conflict[1] && t.genus[1] == "Alpha" && t.genus_agreement[1] == 75
    halves = [ins("Alpha", "Alpha alba"), ins("Alpha", "Alpha alba"),
             ins("Beta", "Alpha alba"), ins("Beta", "Alpha alba")]
    p = bu(fill(99.0, 4), halves)
    @test length(p) == 1 && only(p).rank == :species
    t = lca_assign(hit_table([(fill(99.0, 4), halves)]); direction = :bottom_up)
    @test t.lineage_conflict[1] && t.species[1] == "Alpha alba"
    @test all(r -> ismissing(t[1, r]), (:phylum, :class, :order, :family, :genus))
    # Supporting hits without a name at a higher rank also count as a conflict.
    nofam = [ins("Alpha", "Alpha alba"),
             lineage6("Arthropoda", "Insecta", "Plecoptera", missing, "Alpha", "Alpha alba")]
    p = bu([99.0, 99.0], nofam)
    @test only(filter(x -> x.rank == :family, p)).agreement == 50
    @test lca_assign(hit_table([([99.0, 99.0], nofam)]); direction = :bottom_up).lineage_conflict[1]

    # Not assigned: no hits, or all hits below the phylum threshold.
    @test isempty(bu(Float64[], typeof(ins("A", "A b"))[]))
    @test isempty(bu([40.0, 45.0], [ins("Alpha", "Alpha alba"), ins("Alpha", "Alpha alba")]))
    t = lca_assign(hit_table([(Float64[], []), ([40.0], [ins("Alpha", "Alpha alba")])]);
                   direction = :bottom_up)
    @test all(ismissing, t.assigned_rank) && all(ismissing, t.decision_rank)
    @test !any(t.lineage_conflict) && t.n_hits == [0, 1]

    # Invalid direction: refused before anything else is looked at.
    for bad in (:sideways, "bottom_up", nothing, 1)
        @test_throws r"direction" lca_path([NaN], [ins("Alpha", "Alpha alba")]; direction = bad)
        @test_throws r"direction" lca_assign(DataFrame(); direction = bad)
        @test_throws ArgumentError lca_path([99.0], [ins("Alpha", "Alpha alba")]; direction = bad)
    end
end

@testset "Bottom-up LCA: columns and bold_assign" begin
    sets = [(fill(99.0, 4), [ins("Alpha", "Alpha alba"), ins("Alpha", "Alpha alba"),
                             ins("Alpha", "Alpha alba"), ins("Beta", "Alpha alba")]),
            ([100.0, 97.5], [ins("Alpha", "Alpha alba"), ins("Alpha", "Alpha beta")])]
    hits = hit_table(sets)
    old = legacy_lca_assign(hits)
    top = lca_assign(hits)
    @test propertynames(top) == propertynames(old) && isequal(top, old)
    @test isequal(lca_assign(hits; direction = :top_down), old)
    upw = lca_assign(hits; direction = :bottom_up)
    @test propertynames(upw) == vcat(propertynames(old), [:decision_rank, :lineage_conflict])
    @test isequal(upw.decision_rank, upw.assigned_rank) && upw.lineage_conflict == [true, false]
    @test eltype(upw.lineage_conflict) == Bool
    asvs = DataFrame(id = 1:3, ASV = [seq_of(1), seq_of(2), missing])
    t = bold_assign(asvs; hits, direction = :bottom_up)
    @test isequal(t.decision_rank, ["species", "genus", missing])
    @test isequal(t.lineage_conflict, [true, false, missing])
    @test propertynames(bold_assign(asvs; hits)) == propertynames(t)[1:end-2]
    @test_throws r"decision_rank" bold_assign(hcat(asvs, DataFrame(decision_rank = 1:3)); hits,
                                              direction = :bottom_up)
    @test hasproperty(bold_assign(hcat(asvs, DataFrame(decision_rank = 1:3)); hits), :taxon)
    # Checked before any request to the identification engine.
    silent = FakeEngine(catalog())
    @test_throws r"direction" bold_assign(DataFrame(ASV = [S1]); engine = fake_engine(silent),
                                          direction = :sideways)
    @test isempty(silent.log)
end

@testset "Default direction unchanged since 63e5d70" begin
    rng = TestRNG(0x5eed)
    n_sets = 2500
    n_nonempty = 0
    for _ in 1:n_sets
        sims, lins = random_hits!(rng)
        opts = random_options!(rng)
        limits = BOLDassignment._check_thresholds(opts.thresholds, BU_RANKS)
        want = legacy_lca_path(sims, lins, opts.min_agreement, limits, opts.top_delta,
                               opts.ignore_missing, is_valid_name, BU_RANKS)
        got = lca_path(sims, lins; opts...)
        got == want || @error "Top-down difference" sims lins opts got want
        @test got == want
        @test lca_path(sims, lins; direction = :top_down, opts...) == want
        n_nonempty += !isempty(want)
    end
    @test n_nonempty > 500                     # the sets exercise real assignments
    # Whole tables, several option sets.
    for _ in 1:8
        sets = [random_hits!(rng) for _ in 1:150]
        opts = random_options!(rng)
        hits = hit_table(sets)
        @test isequal(lca_assign(hits; opts...), legacy_lca_assign(hits; opts...))
    end
end

@testset "Bottom-up LCA against an independent oracle" begin
    rng = TestRNG(0xb0770)
    n_sets = 2500
    n_conflicts = 0
    n_assigned = 0
    for _ in 1:n_sets
        sims, lins = random_hits!(rng)
        opts = random_options!(rng)
        want, conflict = oracle_bottom_up(sims, lins, opts)
        got = lca_path(sims, lins; direction = :bottom_up, opts...)
        got == want || @error "Bottom-up difference" sims lins opts got want
        @test got == want
        _, c = BOLDassignment._lca_path_bottom_up(sims, lins, opts.min_agreement,
                                                  BOLDassignment._check_thresholds(opts.thresholds, BU_RANKS),
                                                  opts.top_delta, opts.ignore_missing,
                                                  is_valid_name, BU_RANKS)
        @test c == conflict
        n_conflicts += conflict
        n_assigned += !isempty(want)
    end
    @test n_assigned > 1000 && n_conflicts > 50
    # Tables: decision rank and conflict flags against the oracle.
    for _ in 1:6
        sets = [random_hits!(rng) for _ in 1:150]
        opts = random_options!(rng)
        t = lca_assign(hit_table(sets); direction = :bottom_up, opts...)
        expected = [oracle_bottom_up(s, l, opts) for (s, l) in sets]
        @test t.lineage_conflict == [c for (_, c) in expected]
        @test isequal(t.decision_rank,
                      [isempty(p) ? missing : String(last(p).rank) for (p, _) in expected])
        @test isequal(t.taxon, [isempty(p) ? missing : last(p).taxon for (p, _) in expected])
    end
end
