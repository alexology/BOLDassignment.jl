# Top-hit selection of BOLDigger3 (Buchner & Leese; MIT licence), translated
# from boldigger3 3.0.2, `select_top_hit.py` (`clean_dataframe`,
# `get_threshold`, `move_threshold_up`, `find_top_hit`, `flag_hits`).

"""
    BOLDIGGER3_THRESHOLDS

Similarity thresholds of BOLDigger3 for species, genus, family, order and
class, plus the virtual 50 % it appends for hits known only to phylum.
"""
const BOLDIGGER3_THRESHOLDS = (97, 95, 90, 85, 75, 50)

const _BD3_LEVEL_OF = ("species", "genus", "family", "order", "class", "phylum")  # by threshold
const _BD3_FLAG1 = r"BOLD|ID|Tree|BIN"
# Python's string.punctuation and string.digits, except '-' (clean_dataframe)
const _BD3_INVALID = Set{Char}("!\"#\$%&'()*+,./:;<=>?@[\\]^_`{|}~0123456789")

"""
    boldigger3_name(x) -> Union{Missing,String}

A taxon name as BOLDigger3's `clean_dataframe` leaves it: `missing` for a
missing, empty or `"None"` value and for any name containing an ASCII digit
or punctuation character other than `-` (e.g. `"Baetis sp. 1"`,
`"Gammarus cf. pulex"`).
"""
function boldigger3_name(x)
    (x === missing || x === nothing) && return missing
    s = string(x)
    (isempty(s) || s == "None") && return missing
    any(in(_BD3_INVALID), s) ? missing : s
end

# clean_dataframe: metadata columns lose None/"None"/""; other columns ""
_bd3_metadata(x) = (x === missing || x === nothing || x == "" || x == "None") ? missing : string(x)
_bd3_text(x) = (x === missing || x === nothing || x == "") ? missing : string(x)

const _BD3_COLUMNS = (:id, RANKS..., :pct_identity, :status, :records, :records_ratio,
                      :selected_level, :BIN, :flags, :tied_groups)

"""
    boldigger3_select(hits::AbstractDataFrame; id_col = :id,
                      thresholds = BOLDIGGER3_THRESHOLDS, ties = :first) -> DataFrame

BOLDigger3's top-hit selection, one row per sequence (in order of first
appearance in `hits`), with BOLDigger3's output columns: `id`, the ranks
`phylum` … `species`, `pct_identity`, `status`, `records`, `records_ratio`,
`selected_level`, `BIN`, `flags`, and `tied_groups` (added here: number of
taxa sharing the highest record count at the selected level; with more than
one, BOLDigger3's own choice depends on numpy's sort).

`hits` has one row per hit, in BOLDigger3's order (by decreasing
`pct_identity` within each sequence), with the columns `id_col`, the ranks of
[`RANKS`](@ref), `pct_identity`, `status` (`"public"`/`"private"`), `bin_uri`
and optionally `identification_method`; names are cleaned with
[`boldigger3_name`](@ref). A sequence whose hits contain `"no-match"` is
reported as in BOLDigger3. [`boldigger3_hits`](@ref) builds this table from
[`bold_identify`](@ref) results.

The rule, as in BOLDigger3: the highest similarity picks the starting level
(`thresholds`: species, genus, family, order, class, phylum); the hits at or
above its threshold are grouped by their names down to that level, groups
without a name at that level are dropped (moving up one level if none is
left) and the most numerous group wins. Its hits (all of them, also those
below the threshold) give the reported record, the BINs (species level only)
and the flags: `1` all identification methods name BOLD, ID, Tree or BIN;
`2` the winner has under 90 % of the grouped records; `3` all its records are
private; `4` a single record; `5` several BINs. `ties = :first` takes the
first of tied groups in hit order.
"""
function boldigger3_select(hits::AbstractDataFrame; id_col = :id,
                           thresholds = BOLDIGGER3_THRESHOLDS, ties::Symbol = :first)
    ties === :first || throw(ArgumentError("ties must be :first"))
    thr = Float64.(collect(thresholds))
    length(thr) == 6 && issorted(thr; rev = true) ||
        throw(ArgumentError("thresholds needs 6 decreasing values (species … phylum)"))
    for c in (id_col, RANKS..., :pct_identity, :status, :bin_uri)
        hasproperty(hits, c) || throw(ArgumentError("hits needs a column `$c`"))
    end
    ids = string.(hits[!, id_col])
    tax = Matrix{Union{Missing,String}}(undef, nrow(hits), length(RANKS))
    for (j, r) in enumerate(RANKS)
        tax[:, j] .= boldigger3_name.(hits[!, r])
    end
    raw = [[hits[!, r] for r in RANKS]; [hits[!, :status], hits[!, :bin_uri]]]
    pct = Float64.(hits[!, :pct_identity])
    status = _bd3_text.(hits[!, :status])
    bins = _bd3_text.(hits[!, :bin_uri])
    method = hasproperty(hits, :identification_method) ? _bd3_metadata.(hits[!, :identification_method]) :
             fill(missing, nrow(hits))
    rows = Dict{String,Vector{Int}}(); order = String[]
    for (i, id) in enumerate(ids)
        haskey(rows, id) || push!(order, id)
        push!(get!(rows, id, Int[]), i)
    end
    out = NamedTuple[]
    for id in order
        r = _bd3_top_hit(id, rows[id], tax, raw, pct, status, bins, method, thr)
        r === nothing || push!(out, r)
    end
    df = DataFrame(; (c => Any[getproperty(r, c) for r in out] for c in _BD3_COLUMNS)...)
    for c in _BD3_COLUMNS
        df[!, c] = identity.(df[!, c])  # narrow the element types
    end
    return df
end

function _bd3_result(id, tax_row, pct, status, records, ratio, level, bin, flags, ties)
    (; id, (RANKS[j] => tax_row[j] for j in eachindex(RANKS))..., pct_identity = pct,
       status = status, records = records, records_ratio = ratio, selected_level = level,
       BIN = bin, flags = flags, tied_groups = ties)
end

function _bd3_top_hit(id, idx, tax, raw, pct, status, bins, method, thr)
    nolevels = fill(missing, length(RANKS))
    # get_threshold: a "no-match" anywhere in the hits of the sequence
    if any(col -> any(i -> isequal(col[i], "no-match"), idx), raw)
        i = findfirst(i -> isequal(tax[i, end], "no-match"), idx)
        i === nothing && return nothing  # BOLDigger3 returns an empty row set
        k = idx[i]
        return _bd3_result(id, tax[k, :], pct[k], status[k], 0, missing, missing, missing, "||||", 0)
    end
    best = maximum(pct[i] for i in idx)
    t = something(findfirst(x -> best >= x, thr[1:5]), 6)
    while true
        li = length(RANKS) + 1 - t          # rank index of the level (species = 6)
        counts = Dict{Vector{Union{Missing,String}},Int}(); keys_order = Vector{Union{Missing,String}}[]
        for i in idx
            pct[i] >= thr[t] || continue
            key = tax[i, 1:li]
            haskey(counts, key) || push!(keys_order, key)
            counts[key] = get(counts, key, 0) + 1
        end
        filter!(k -> k[li] !== missing, keys_order)  # dropna at the selected level
        if isempty(keys_order)
            if t == 6
                k = idx[1]
                return _bd3_result(id, nolevels, pct[k], status[k], 0, 0.0, missing, missing, "||||", 0)
            end
            t += 1
            continue
        end
        total = sum(counts[k] for k in keys_order)
        top = maximum(counts[k] for k in keys_order)
        tied = count(k -> counts[k] == top, keys_order)
        winner = keys_order[findfirst(k -> counts[k] == top, keys_order)]
        ratio = top / total
        # query on the named levels of the winner, over all hits of the sequence
        top_hits = [i for i in idx if all(j -> winner[j] === missing ||
                                               (tax[i, j] !== missing && tax[i, j] == winner[j]), 1:li)]
        bin = t == 1 ? join(unique(skipmissing(bins[i] for i in top_hits)), "|") : ""
        k = top_hits[1]
        row = Vector{Union{Missing,String}}(tax[k, :])
        t == 1 || (row[li+1:end] .= missing)
        methods = collect(skipmissing(method[i] for i in top_hits))
        flags = ["", "", "", "", ""]
        !isempty(methods) && all(m -> occursin(_BD3_FLAG1, m), methods) && (flags[1] = "1")
        ratio < 0.9 && (flags[2] = "2")
        all(i -> isequal(status[i], "private"), top_hits) && (flags[3] = "3")
        length(top_hits) == 1 && (flags[4] = "4")
        length(split(bin, '|')) > 1 && (flags[5] = "5")
        return _bd3_result(id, row, pct[k], status[k], top, ratio, _BD3_LEVEL_OF[t], bin,
                           join(flags, '|'), tied)
    end
end

"""
    boldigger3_hits(hits::AbstractDataFrame; ids = nothing, metadata = nothing) -> DataFrame

The input of [`boldigger3_select`](@ref) from a [`bold_identify`](@ref) table,
as BOLDigger3 builds it: one row per hit with `id` (the sequence, or
`ids[sequence]` when `ids` maps sequences to names), the ranks,
`pct_identity` (= `similarity`), `process_id`, `bin_uri`, hits sorted by
decreasing similarity within each sequence (stable: BOLD's order on ties) and
a `"no-match"` row for a sequence without hits.

Hits with equal similarity keep BOLD's order; BOLDigger3 receives them in
the order of its DuckDB query, which is not defined for equal values. This
only matters for the order of the BINs in `BIN` and for tied groups
(`tied_groups > 1`).

`metadata` (optional) plays the role of BOLDigger3's local database: a table
with `processid` and `identification_method`, one row per record and marker.
Each hit is joined to every row of its process ID (so a record with several
markers counts several times, as in BOLDigger3), `status` is `"public"` when
a row was found and `"private"` otherwise. Without `metadata`, `status` is
BOLD's `record_status` and `identification_method` is missing, so flags 1
and 3 may differ from BOLDigger3.
"""
function boldigger3_hits(hits::AbstractDataFrame; ids = nothing, metadata = nothing)
    name(s) = ids === nothing ? s : string(ids[s])
    meta = Dict{String,Vector{Union{Missing,String}}}()
    if metadata !== nothing
        for (p, m) in zip(metadata[!, :processid], metadata[!, :identification_method])
            p === missing && continue
            push!(get!(meta, string(p), Union{Missing,String}[]), m === missing ? missing : string(m))
        end
    end
    out = (id = String[], process_id = Union{Missing,String}[], bin_uri = Union{Missing,String}[],
           pct_identity = Float64[], status = Union{Missing,String}[],
           identification_method = Union{Missing,String}[])
    taxa = Dict(r => Union{Missing,String}[] for r in RANKS)
    rowsof = Dict{String,Vector{Int}}(); order = String[]
    for (i, s) in enumerate(hits[!, :sequence])
        haskey(rowsof, s) || push!(order, s)
        push!(get!(rowsof, s, Int[]), i)
    end
    for s in order
        idx = rowsof[s]
        real = [i for i in idx if hits[i, :process_id] !== missing]
        if isempty(real)
            push!(out.id, name(s)); push!(out.process_id, missing); push!(out.bin_uri, missing)
            # BOLDigger3 joins the no-match row too: no process ID, hence "private"
            push!(out.pct_identity, 0.0); push!(out.status, metadata === nothing ? missing : "private")
            push!(out.identification_method, missing)
            foreach(r -> push!(taxa[r], "no-match"), RANKS)
            continue
        end
        sort!(real; by = i -> -hits[i, :similarity])  # stable
        for i in real
            pid = string(hits[i, :process_id])
            joined = metadata === nothing ? [missing] : get(meta, pid, nothing)
            st = metadata === nothing ? hits[i, :record_status] : (joined === nothing ? "private" : "public")
            for m in something(joined, [missing])
                push!(out.id, name(s)); push!(out.process_id, pid); push!(out.bin_uri, hits[i, :bin_uri])
                push!(out.pct_identity, hits[i, :similarity]); push!(out.status, st)
                push!(out.identification_method, m)
                foreach(r -> push!(taxa[r], hits[i, r]), RANKS)
            end
        end
    end
    df = DataFrame(out)
    for r in RANKS
        df[!, r] = taxa[r]
    end
    return select!(df, :id, RANKS..., :pct_identity, :process_id, :bin_uri, :status, :identification_method)
end
