mutable struct _Chunk
    seqs::Vector{String}
    file::Union{Nothing,String}
    sub_id::Union{Nothing,String}
    submitted_at::Float64
    checked_at::Float64
    submissions::Int
end
_Chunk(seqs, file) = _Chunk(seqs, file, nothing, -Inf, -Inf, 0)

"""
    bold_identify(sequences; db = :animal_public, mode = :rapid, engine = IDEngine(),
                  chunk_size = nothing, checkpoint = nothing) -> DataFrame

Identify `sequences` with the BOLD identification engine and return all hits
as a long table.

Sequences are normalized with [`normalize_sequence`](@ref), submitted once
each and grouped in chunks of at most `chunk_size` (default: the maximum
documented for `mode`, see [`BOLD_MODES`](@ref)) after sorting, so chunks do
not depend on the input order. `db` is a key of [`BOLD_DATABASES`](@ref) or
an API code. With `checkpoint` set to a directory, the raw response of each
chunk is saved there and reused by later calls with the same sequences,
database and mode, so an interrupted run resumes where it stopped. A
checkpoint file that cannot be read back is reported and its chunk submitted
again.

Columns: `sequence` (normalized), `process_id`, `bin_uri`, `record_status`,
`similarity` (percent, `100 - pdist`) and one column per rank in
[`RANKS`](@ref). A sequence without hits has one row with `missing` hit
fields. Rows follow the sorted order of the sequences.
"""
function bold_identify(sequences::AbstractVector; db = :animal_public, mode::Symbol = :rapid,
                       engine::IDEngine = IDEngine(), chunk_size::Union{Nothing,Integer} = nothing,
                       checkpoint::Union{Nothing,AbstractString} = nothing)
    code = _database_code(db)
    params = _mode(mode)
    size = chunk_size === nothing ? params.max_sequences : Int(chunk_size)
    1 <= size <= params.max_sequences ||
        throw(ArgumentError("chunk_size must be between 1 and $(params.max_sequences) for mode :$mode"))
    unique_seqs = sort!(unique!(String[normalize_sequence(s) for s in sequences if !ismissing(s)]))
    if checkpoint !== nothing
        ispath(checkpoint) && !isdir(checkpoint) &&
            throw(ArgumentError("checkpoint must be a directory: $checkpoint"))
        mkpath(checkpoint)
    end

    chunks = [_Chunk(unique_seqs[r], nothing) for r in Iterators.partition(eachindex(unique_seqs), size)]
    tables = Vector{Union{Nothing,DataFrame}}(nothing, length(chunks))
    pending = Int[]
    for (i, c) in enumerate(chunks)
        if checkpoint !== nothing
            c.file = joinpath(checkpoint, _chunk_key(code, params, c.seqs) * ".jsonl")
            if isfile(c.file)
                tables[i] = try
                    _hits_table(c.seqs, read(c.file))
                catch err
                    err isa InterruptException && rethrow()
                    @warn "Unreadable checkpoint $(c.file); submitting its chunk again" exception = err
                    nothing
                end
                tables[i] === nothing || continue
            end
        end
        push!(pending, i)
    end
    reused = length(chunks) - length(pending)
    reused > 0 && _log(engine, "Reusing $reused checkpointed chunk(s)")

    _run_chunks!(tables, chunks, pending, engine, code, params)

    isempty(tables) && return _hits_table(String[], UInt8[])
    return reduce(vcat, tables)
end

function _chunk_key(code, params, seqs)
    ctx = SHA.SHA256_CTX()
    SHA.update!(ctx, codeunits(string(code, '|', params.mi, '|', params.mo, '|', params.maxh, '\n')))
    for s in seqs
        SHA.update!(ctx, codeunits(s))
        SHA.update!(ctx, UInt8['\n'])
    end
    return bytes2hex(SHA.digest!(ctx))
end

function _run_chunks!(tables, chunks, pending, engine, code, params)
    active = Int[]
    queue = copy(pending)
    total = length(pending)
    done = 0
    while !isempty(queue) || !isempty(active)
        if !isempty(queue) && length(active) < engine.max_active
            i = popfirst!(queue)
            c = chunks[i]
            c.submissions > engine.max_resubmits &&
                throw(ErrorException("chunk $i was submitted $(c.submissions) times without results"))
            ids = string.(eachindex(c.seqs))
            c.sub_id = _submit(engine, code, params, _fasta_bytes(c.seqs, ids))
            c.submissions += 1
            c.submitted_at = engine.time()
            c.checked_at = c.submitted_at
            push!(active, i)
            _log(engine, "Submitted chunk $i ($(length(c.seqs)) sequences): $(c.sub_id)")
            continue
        end
        next_check = minimum(chunks[i].checked_at + engine.status_wait for i in active)
        _wait_until(engine, next_check)
        for i in copy(active)
            c = chunks[i]
            engine.time() < c.checked_at + engine.status_wait && continue
            result = _poll(engine, c)
            c.checked_at = engine.time()
            if result !== nothing
                body, table = result
                c.file === nothing || _save_atomic(c.file, body)
                tables[i] = table
                filter!(!=(i), active)
                done += 1
                _log(engine, "Chunk $i finished ($done/$total)")
            elseif c.checked_at - c.submitted_at > engine.timeout
                _log(engine, "Chunk $i timed out; resubmitting")
                filter!(!=(i), active)
                push!(queue, i)
            end
        end
    end
    return tables
end

# `(body, table)` once the submission is complete and its results parse with
# every sequence present; `nothing` otherwise.
function _poll(engine, c::_Chunk)
    _finished(engine, c.sub_id) || return nothing
    body = _results(engine, c.sub_id)
    body === nothing && return nothing
    table = try
        _hits_table(c.seqs, body)
    catch err
        err isa InterruptException && rethrow()
        _log(engine, "Results of $(c.sub_id) not usable yet: $(sprint(showerror, err))")
        return nothing
    end
    return body, table
end

# Write to a unique temporary file in the same directory, then rename over
# the destination; the temporary file is removed on failure.
function _save_atomic(file::AbstractString, body::Vector{UInt8})
    tmp = tempname(dirname(file); cleanup = false) * ".partial"
    try
        write(tmp, body)
        Base.Filesystem.rename(tmp, file)
    catch
        rm(tmp; force = true)
        rethrow()
    end
    return nothing
end

_text(x) = (x === nothing || ismissing(x)) ? missing :
           (s = strip(string(x)); isempty(s) || s == "None" ? missing : String(s))

function _similarity(rec, key)
    pdist = rec isa JSON3.Object ? get(rec, :pdist, nothing) : nothing
    pdist isa Real && isfinite(pdist) && 0 <= pdist <= 100 ||
        throw(ErrorException("invalid pdist $(repr(pdist)) for hit $key"))
    return 100.0 - Float64(pdist)
end

# Hit table of one chunk from the JSON Lines body; seqids are positions in `seqs`.
function _hits_table(seqs::AbstractVector{String}, body::Vector{UInt8})
    cols = (sequence = String[], process_id = Union{Missing,String}[],
            bin_uri = Union{Missing,String}[], record_status = Union{Missing,String}[],
            similarity = Union{Missing,Float64}[])
    taxa = Dict(r => Union{Missing,String}[] for r in RANKS)
    function push_row!(seq, pid, bin, status, sim, tax)
        push!(cols.sequence, seq); push!(cols.process_id, pid); push!(cols.bin_uri, bin)
        push!(cols.record_status, status); push!(cols.similarity, sim)
        for r in RANKS
            push!(taxa[r], tax isa JSON3.Object ? _text(get(tax, r, nothing)) : missing)
        end
    end
    seen = falses(length(seqs))
    for line in eachline(IOBuffer(body))
        isempty(strip(line)) && continue
        obj = JSON3.read(line)
        k = tryparse(Int, string(get(obj, :seqid, "")))
        (k !== nothing && 1 <= k <= length(seqs)) ||
            throw(ErrorException("unexpected seqid $(repr(get(obj, :seqid, nothing))) in BOLD results"))
        seen[k] && throw(ErrorException("seqid $k repeated in BOLD results"))
        seen[k] = true
        results = get(obj, :results, nothing)
        if results === nothing || isempty(results)
            push_row!(seqs[k], missing, missing, missing, missing, nothing)
            continue
        end
        results isa JSON3.Object || throw(ErrorException("unexpected results for seqid $k"))
        for (key, rec) in pairs(results)
            fields = split(string(key), '|')
            pid = _text(fields[1])
            pid === missing && throw(ErrorException("hit without process ID for seqid $k"))
            push_row!(seqs[k], pid,
                      length(fields) >= 3 ? _text(fields[3]) : missing,
                      length(fields) >= 5 ? _text(fields[5]) : missing,
                      _similarity(rec, key), get(rec, :taxonomy, nothing))
        end
    end
    all(seen) || throw(ErrorException("BOLD results miss $(count(!, seen)) sequence(s)"))
    df = DataFrame(cols)
    for r in RANKS
        df[!, r] = taxa[r]
    end
    return df
end
