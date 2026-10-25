mutable struct _Chunk
    index::Vector{Int}   # positions in the sorted unique sequences
    seqs::Vector{String}
    sub_id::Union{Nothing,String}
    submitted_at::Float64
    checked_at::Float64
    submissions::Int
end
_Chunk(index, seqs) = _Chunk(index, seqs, nothing, -Inf, -Inf, 0)

"""
    bold_identify(sequences; db = :animal_public, mode = :rapid, engine = IDEngine(),
                  chunk_size = nothing, checkpoint = nothing) -> DataFrame

Identify `sequences` with the BOLD identification engine and return all hits
as a long table.

Sequences are normalized with [`normalize_sequence`](@ref), submitted once
each and grouped in chunks of at most `chunk_size` (default: the maximum
documented for `mode`, see [`BOLD_MODES`](@ref)) after sorting. With
`checkpoint` set to a directory, the raw result of each sequence is saved
there (one file per sequence, database and mode) and reused by later calls,
so an interrupted run resumes where it stopped and adding sequences submits
only the new ones. A checkpoint file that cannot be read back is reported
and its sequence submitted again. `db` is a key of [`BOLD_DATABASES`](@ref)
or an API code.

If a finished submission lacks the result of some sequences, they are
reported with a warning, returned without hits and not saved in the
checkpoint, so a later call submits them again.

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

    lines = Vector{Union{Nothing,String}}(nothing, length(unique_seqs))
    files = checkpoint === nothing ? nothing :
            [joinpath(checkpoint, _sequence_key(code, params, s) * ".jsonl") for s in unique_seqs]
    pending = Int[]
    for (i, s) in enumerate(unique_seqs)
        if files !== nothing && isfile(files[i])
            line = try
                l = read(files[i], String)
                _hits_table([s], [l])  # validates the saved result
                l
            catch err
                err isa InterruptException && rethrow()
                @warn "Unreadable checkpoint $(files[i]); submitting its sequence again" exception = err
                nothing
            end
            line === nothing || (lines[i] = line; continue)
        end
        push!(pending, i)
    end
    reused = length(unique_seqs) - length(pending)
    reused > 0 && _log(engine, "Reusing $reused checkpointed sequence(s)")

    chunks = [_Chunk(pending[r], unique_seqs[pending[r]]) for r in Iterators.partition(eachindex(pending), size)]
    _run_chunks!(lines, files, chunks, engine, code, params, checkpoint)
    return _hits_table(unique_seqs, lines)
end

function _sequence_key(code, params, seq)
    ctx = SHA.SHA256_CTX()
    SHA.update!(ctx, codeunits(string(code, '|', params.mi, '|', params.mo, '|', params.maxh, '\n')))
    SHA.update!(ctx, codeunits(seq))
    return bytes2hex(SHA.digest!(ctx))
end

function _run_chunks!(lines, files, chunks, engine, code, params, checkpoint)
    active = Int[]
    queue = collect(eachindex(chunks))
    total = length(chunks)
    done = 0
    while !isempty(queue) || !isempty(active)
        if !isempty(queue) && length(active) < engine.max_active
            i = popfirst!(queue)
            c = chunks[i]
            if c.submissions > engine.max_resubmits
                msg = "chunk $i ($(length(c.seqs)) sequences) was submitted $(c.submissions) times " *
                      "without usable results; last submission $(c.sub_id)."
                checkpoint === nothing ||
                    (msg *= " Finished sequences are saved in $checkpoint; call again to resume.")
                throw(ErrorException(msg))
            end
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
                for (k, line) in zip(c.index, result)
                    lines[k] = line
                    line === nothing || files === nothing || _save_atomic(files[k], line)
                end
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
    return lines
end

# Result line of each sequence of the chunk (`nothing` for sequences missing
# from a finished submission) once the results are complete and valid;
# `nothing` while they are not available yet.
function _poll(engine, c::_Chunk)
    _finished(engine, c.sub_id) || return nothing
    body = _results(engine, c.sub_id)
    body === nothing && return nothing
    lines = try
        l = _result_lines(c.seqs, body; allow_missing = true)
        _hits_table(c.seqs, l)  # validates every hit before anything is saved
        l
    catch err
        err isa InterruptException && rethrow()
        _log(engine, "Results of $(c.sub_id) not usable yet: $(sprint(showerror, err))")
        return nothing
    end
    n = count(isnothing, lines)
    n > 0 && @warn "BOLD submission $(c.sub_id) finished without the result of $n sequence(s); " *
                   "they are reported without hits and not checkpointed"
    return lines
end

# Write to a unique temporary file in the same directory, then rename over
# the destination; the temporary file is removed on failure.
function _save_atomic(file::AbstractString, body)
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

# Splits a JSON Lines body into the line of each sequence; seqids are positions
# in `seqs`. Missing sequences are an error unless `allow_missing`.
function _result_lines(seqs::AbstractVector{String}, body::Vector{UInt8}; allow_missing::Bool = false)
    lines = Vector{Union{Nothing,String}}(nothing, length(seqs))
    for line in eachline(IOBuffer(body))
        isempty(strip(line)) && continue
        obj = JSON3.read(line)
        k = tryparse(Int, string(get(obj, :seqid, "")))
        (k !== nothing && 1 <= k <= length(seqs)) ||
            throw(ErrorException("unexpected seqid $(repr(get(obj, :seqid, nothing))) in BOLD results"))
        lines[k] === nothing || throw(ErrorException("seqid $k repeated in BOLD results"))
        lines[k] = line
    end
    n = count(isnothing, lines)
    n > 0 && !allow_missing && throw(ErrorException("BOLD results miss $n sequence(s)"))
    return lines
end

_hits_table(seqs::AbstractVector{String}, body::Vector{UInt8}) = _hits_table(seqs, _result_lines(seqs, body))

# Hit table from the result line of each sequence (`nothing`: no hits). The
# seqid inside a line is not used here.
function _hits_table(seqs::AbstractVector{String}, lines::AbstractVector)
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
    for (k, line) in enumerate(lines)
        results = line === nothing ? nothing : get(JSON3.read(line), :results, nothing)
        if results === nothing || isempty(results)
            push_row!(seqs[k], missing, missing, missing, missing, nothing)
            continue
        end
        results isa JSON3.Object || throw(ErrorException("unexpected results for sequence $k"))
        for (key, rec) in pairs(results)
            fields = split(string(key), '|')
            pid = _text(fields[1])
            pid === missing && throw(ErrorException("hit without process ID for sequence $k"))
            push_row!(seqs[k], pid,
                      length(fields) >= 3 ? _text(fields[3]) : missing,
                      length(fields) >= 5 ? _text(fields[5]) : missing,
                      _similarity(rec, key), get(rec, :taxonomy, nothing))
        end
    end
    df = DataFrame(cols)
    for r in RANKS
        df[!, r] = taxa[r]
    end
    return df
end
