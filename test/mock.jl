# Fake identification engine, in the response format observed on 2026-10-23.
# Answers from `catalog` (sequence => vector of (process_id, bin, pdist,
# taxonomy::Dict)) and records every request.
mutable struct FakeEngine
    catalog::Dict{String,Vector{Any}}
    submissions::Dict{String,Vector{Pair{String,String}}}
    polls_before_done::Int
    polls::Dict{String,Int}
    log::Vector{Tuple{Float64,String,String}}
    clock::Float64
    fail_next::Vector{Int}
end
FakeEngine(catalog; polls_before_done = 1) =
    FakeEngine(catalog, Dict(), polls_before_done, Dict(), [], 0.0, Int[])

function fake_engine(fake::FakeEngine; kwargs...)
    IDEngine(; base_url = "https://id.test", request = (m, u, p) -> fake_request(fake, m, u, p),
             sleep = dt -> (fake.clock += dt), time = () -> fake.clock, verbose = false,
             throttle = BOLDassignment._Throttle(), kwargs...)
end

function fake_request(fake::FakeEngine, method, url, payload)
    push!(fake.log, (fake.clock, method, url))
    fake.clock += 0.1
    if !isempty(fake.fail_next)
        return (popfirst!(fake.fail_next), UInt8[])
    end
    if method == "POST"
        records = Pair{String,String}[]
        for rec in FASTA.Reader(IOBuffer(payload))
            push!(records, String(FASTA.identifier(rec)) => String(FASTA.sequence(rec)))
        end
        sub_id = "sub$(length(fake.submissions) + 1)"
        fake.submissions[sub_id] = records
        return (200, Vector{UInt8}(codeunits("{\"sub_id\":\"$sub_id\",\"parameters\":{}}")))
    end
    sub_id = last(split(url, '/'))
    haskey(fake.submissions, sub_id) || return (404, UInt8[])
    n = length(fake.submissions[sub_id])
    if occursin("/status/", url)
        fake.polls[sub_id] = get(fake.polls, sub_id, 0) + 1
        done = fake.polls[sub_id] > fake.polls_before_done
        body = done ? "{\"status\":\"complete\"}" :
                      "{\"queued\":0,\"processing\":$n,\"completed\":0}"
        return (200, Vector{UInt8}(codeunits(body)))
    end
    io = IOBuffer()
    for (id, seq) in fake.submissions[sub_id]
        results = Dict{String,Any}()
        for (pid, bin, pdist, tax) in get(fake.catalog, seq, [])
            results["$pid|COI-5P|$bin|123456|public|0"] =
                Dict("pdist" => pdist, "pident" => 100 - pdist,
                     "taxonomy" => merge(tax, Dict("taxid_count" => "1")))
        end
        println(io, JSON3.write(Dict("seqid" => id, "results" => results)))
    end
    return (200, take!(io))
end

lineage(phylum, class, order, family, genus, species) =
    Dict("phylum" => phylum, "class" => class, "order" => order, "family" => family,
         "genus" => genus, "species" => species)
beetle(genus, species) = lineage("Arthropoda", "Insecta", "Coleoptera", "Scarabaeidae", genus, species)
