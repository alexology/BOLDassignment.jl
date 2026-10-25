const _PORTAL_URL = "https://portal.boldsystems.org"
# The portal `query` parameter accepts at most 250 characters.
const _PORTAL_QUERY_MAX = 250

"""
    bold_records(process_ids; fields = nothing, engine = IDEngine(), base_url =
                 "https://portal.boldsystems.org") -> DataFrame

Retrieve public BOLD records (BCDM fields) for `process_ids` through the
portal API, e.g. to add collection data to the hits of
[`bold_identify`](@ref). Only the columns in `fields` are kept when given;
nested values (such as `coord`) are returned as JSON text. Unknown or private
process IDs are absent from the result. Requests use the spacing, retries and
transport of `engine`.
"""
function bold_records(process_ids::AbstractVector; fields = nothing, engine::IDEngine = IDEngine(),
                      base_url::AbstractString = _PORTAL_URL)
    ids = unique!(String[strip(string(p)) for p in process_ids if !ismissing(p)])
    filter!(!isempty, ids)
    records = Dict{Symbol,Any}[]
    for batch in _portal_batches(ids)
        query = join(("ids:processid:" * p for p in batch), ';')
        _, body = _call(engine, "GET",
                        string(base_url, "/api/query?query=", HTTP.escapeuri(query), "&extent=full"))
        qid = String(JSON3.read(body)[:query_id])
        start = 0
        while true
            _, body = _call(engine, "GET",
                            string(base_url, "/api/documents/", HTTP.escapeuri(qid),
                                   "?length=100&start=", start))
            data = JSON3.read(body)[:data]
            for rec in data
                push!(records, Dict{Symbol,Any}(k => _portal_value(v) for (k, v) in pairs(rec)))
            end
            length(data) < 100 && break
            start += 100
        end
    end
    columns = fields === nothing ? _union_keys(records) : Symbol.(collect(fields))
    df = DataFrame()
    for c in columns
        # Narrowest element type (e.g. Union{Missing,String}) instead of Any.
        df[!, c] = identity.([get(r, c, missing) for r in records])
    end
    return df
end

function _portal_batches(ids)
    batches = Vector{String}[]
    current = String[]
    len = 0
    for p in ids
        term = length("ids:processid:") + length(p)
        term > _PORTAL_QUERY_MAX && throw(ArgumentError("process ID too long: $p"))
        # One more character for the ';' separator when the batch is not empty.
        if !isempty(current) && len + 1 + term > _PORTAL_QUERY_MAX
            push!(batches, current)
            current = String[]
            len = 0
        end
        len += term + (isempty(current) ? 0 : 1)
        push!(current, p)
    end
    isempty(current) || push!(batches, current)
    return batches
end

_portal_value(::Nothing) = missing
_portal_value(v::Union{AbstractString,Number,Bool}) = v isa AbstractString ? String(v) : v
_portal_value(v) = JSON3.write(v)

function _union_keys(records)
    seen = Set{Symbol}()
    keys_ordered = Symbol[]
    for r in records, k in keys(r)
        k in seen || (push!(seen, k); push!(keys_ordered, k))
    end
    return keys_ordered
end
