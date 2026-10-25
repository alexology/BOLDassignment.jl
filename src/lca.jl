# Epithets and names used by open nomenclature or as placeholders.
const _PLACEHOLDER_EPITHETS = Set(["sp", "spp", "cf", "aff", "nr", "near", "gr", "group",
                                   "complex", "indet", "undet", "nov", "n", "x"])
const _PLACEHOLDER_NAMES = Set(["unknown", "unclassified", "undetermined", "unidentified",
                                "incertae", "none", "na", "null"])
const _SPECIES_NAME = r"^[A-Z][a-z]+(-[a-z]+)?( \([A-Z][a-z]+\))? [a-z]+(-[a-z]+)?( [a-z]+(-[a-z]+)?)?$"
const _TAXON_NAME = r"^[A-Z][a-z]+(-[a-z]+)?$"

"""
    is_valid_name(rank, name) -> Bool

Whether `name` is a usable identification at `rank`. A species must be a
Latin binomial (optionally with subgenus in parentheses or a subspecies
epithet) without open nomenclature: names such as `Onthophagus sp.`,
`Onthophagus sp. 1`, `Onthophagus cf. taurus`, `Onthophagus aff. taurus` or
interim codes (`Diptera sp. BIOUG1234-A01`) are rejected. Other ranks must be
a single capitalized word without digits, and not a placeholder such as
`Unknown`. `missing` is never valid.
"""
function is_valid_name(rank::Symbol, name)
    ismissing(name) && return false
    s = strip(string(name))
    if rank === :species
        occursin(_SPECIES_NAME, s) || return false
        words = split(replace(s, r" \([A-Z][a-z]+\)" => ""))
        return !any(w -> w in _PLACEHOLDER_EPITHETS, words[2:end])
    end
    return occursin(_TAXON_NAME, s) && !(lowercase(s) in _PLACEHOLDER_NAMES)
end

const _ASSIGNMENT_COLUMNS = (:sequence, :n_hits, :best_similarity, :best_process_id,
                             :assigned_rank, :taxon)

function _check_ranks(ranks)
    rs = Symbol[Symbol(r) for r in ranks]
    isempty(rs) && throw(ArgumentError("ranks must not be empty"))
    allunique(rs) || throw(ArgumentError("ranks must be unique"))
    names = vcat(rs, [Symbol(r, :_agreement) for r in rs])
    for r in rs
        r in _ASSIGNMENT_COLUMNS && throw(ArgumentError("rank :$r clashes with an output column"))
    end
    allunique(names) || throw(ArgumentError("rank names clash with the `<rank>_agreement` columns"))
    return Tuple(rs)
end

# Thresholds given by the user override the defaults rank by rank.
function _check_thresholds(thresholds, ranks)
    limits = Dict{Symbol,Float64}(r => DEFAULT_THRESHOLDS[r] for r in ranks
                                  if haskey(DEFAULT_THRESHOLDS, r))
    for (k, v) in pairs(thresholds)
        r = Symbol(k)
        r in ranks || throw(ArgumentError("threshold for :$r, which is not in ranks"))
        (v isa Real && isfinite(v) && 0 <= v <= 100) ||
            throw(ArgumentError("threshold for :$r must be within [0, 100], got $v"))
        limits[r] = v
    end
    for r in ranks
        haskey(limits, r) || throw(ArgumentError("no threshold for rank :$r"))
    end
    return limits
end

function _check_lca_options(min_agreement, top_delta)
    0 < min_agreement <= 100 || throw(ArgumentError("min_agreement must be in (0, 100]"))
    (top_delta >= 0 && !isnan(top_delta)) ||
        throw(ArgumentError("top_delta must be non-negative"))
    return nothing
end

"""
    lca_path(similarity, lineages; min_agreement = 100, thresholds = DEFAULT_THRESHOLDS,
             top_delta = Inf, ignore_missing = true, valid_name = is_valid_name,
             ranks = RANKS)

Lowest common ancestor of the hits of one sequence with a minimum agreement.

`similarity` holds the percent similarity of each hit and `lineages` one
indexable lineage per hit (`lineage[r]` for each rank `r` in `ranks`, `missing`
when unknown). Ranks are visited from the most general to the most specific.
At rank `r` the pool contains the hits with similarity at least
`thresholds[r]` and at least `maximum(similarity) - top_delta`. The most
frequent taxon at `r` among the pool hits that agree with the ranks already
assigned is accepted when it is supported by at least `min_agreement` percent
of the pool; otherwise the assignment stops at the previous rank. Ties for
the most frequent taxon stop the assignment.

With `min_agreement = 100` (default) this is the strict LCA of the pool.
The default `top_delta = Inf` keeps every hit above the rank threshold, the
conservative choice; a finite value restricts the pool to the closest hits.
`thresholds` may list only some ranks: the others keep
[`DEFAULT_THRESHOLDS`](@ref). A rank without a default threshold must be
given one.

Hits not identified at `r`, i.e. with a `missing` name or one rejected by
`valid_name(r, name)` (default [`is_valid_name`](@ref)), leave the pool at
that rank and at the deeper ones: a record of `Onthophagus sp.` or
`Onthophagus cf. taurus` above the species threshold does not take part in
the species decision but still supports the genus. With
`ignore_missing = false` such hits stay in the pool and count against the
agreement.

Returns a vector of `(rank, taxon, agreement, n_support, n_pool)` named tuples,
one per assigned rank (empty when nothing can be assigned).
"""
function lca_path(similarity::AbstractVector{<:Real}, lineages::AbstractVector;
                  min_agreement::Real = 100, thresholds = DEFAULT_THRESHOLDS,
                  top_delta::Real = Inf, ignore_missing::Bool = true,
                  valid_name = is_valid_name, ranks = RANKS)
    ranks = _check_ranks(ranks)
    limits = _check_thresholds(thresholds, ranks)
    _check_lca_options(min_agreement, top_delta)
    length(similarity) == length(lineages) ||
        throw(ArgumentError("similarity and lineages must have the same length"))
    for s in similarity
        (isfinite(s) && 0 <= s <= 100) ||
            throw(ArgumentError("similarity must be finite and within [0, 100], got $s"))
    end
    return _lca_path(similarity, lineages, min_agreement, limits, top_delta, ignore_missing,
                     valid_name, ranks)
end

function _lca_path(similarity, lineages, min_agreement, limits, top_delta, ignore_missing,
                   valid_name, ranks)
    path = NamedTuple{(:rank, :taxon, :agreement, :n_support, :n_pool),
                      Tuple{Symbol,String,Float64,Int,Int}}[]
    isempty(similarity) && return path
    floor_sim = maximum(similarity) - top_delta
    # Hits that agree with every rank assigned so far.
    consistent = trues(length(similarity))
    # Hits not identified at a previous rank, left out of the deeper ranks.
    unidentified = falses(length(similarity))
    counts = Dict{String,Int}()
    for r in ranks
        t = limits[r]
        empty!(counts)
        n_pool = 0
        for i in eachindex(similarity)
            s = similarity[i]
            (s >= t && s >= floor_sim) || continue
            unidentified[i] && continue
            name = lineages[i][r]
            if !valid_name(r, name)
                if ignore_missing
                    unidentified[i] = true
                else
                    n_pool += 1
                end
                continue
            end
            n_pool += 1
            consistent[i] && (counts[name] = get(counts, name, 0) + 1)
        end
        (n_pool == 0 || isempty(counts)) && break
        best = maximum(values(counts))
        count(==(best), values(counts)) == 1 || break
        taxon = first(k for (k, v) in counts if v == best)
        agreement = 100 * best / n_pool
        # Relative tolerance so that e.g. 2/3 passes min_agreement = 66.66666666666667.
        agreement >= min_agreement * (1 - 1e-12) || break
        push!(path, (rank = r, taxon = taxon, agreement = agreement, n_support = best,
                     n_pool = n_pool))
        for i in eachindex(similarity)
            consistent[i] && (consistent[i] = isequal(lineages[i][r], taxon))
        end
    end
    return path
end

"""
    lca_assign(hits::DataFrame; min_agreement = 100, thresholds = DEFAULT_THRESHOLDS,
               top_delta = Inf, ignore_missing = true, valid_name = is_valid_name,
               ranks = RANKS) -> DataFrame

Apply [`lca_path`](@ref) to every sequence of a hit table returned by
[`bold_identify`](@ref).

One row per sequence, in order of first appearance, with columns `sequence`,
`n_hits`, `best_similarity`, `best_process_id`, `assigned_rank`, `taxon`
(deepest assigned taxon), one column per rank with the assigned name and one
`<rank>_agreement` column with the percent of pool hits supporting it.
Unassigned ranks are `missing`; a sequence without hits has `n_hits = 0`.
"""
function lca_assign(hits::AbstractDataFrame; min_agreement::Real = 100,
                    thresholds = DEFAULT_THRESHOLDS, top_delta::Real = Inf,
                    ignore_missing::Bool = true, valid_name = is_valid_name,
                    ranks = RANKS)
    ranks = _check_ranks(ranks)
    limits = _check_thresholds(thresholds, ranks)
    _check_lca_options(min_agreement, top_delta)
    for col in (:sequence, :similarity, :process_id, ranks...)
        hasproperty(hits, col) || throw(ArgumentError("hit table lacks column :$col"))
    end
    out = DataFrame(sequence = String[], n_hits = Int[],
                    best_similarity = Union{Missing,Float64}[],
                    best_process_id = Union{Missing,String}[],
                    assigned_rank = Union{Missing,String}[], taxon = Union{Missing,String}[])
    for r in ranks
        out[!, r] = Union{Missing,String}[]
    end
    for r in ranks
        out[!, Symbol(r, :_agreement)] = Union{Missing,Float64}[]
    end
    groups, order, repeated = _hit_groups(hits, ranks)
    repeated > 0 && @warn "Ignored $repeated repeated hit row(s) (same sequence, process ID and data)"
    for s in order
        idx = groups[s]
        sims = Float64[hits.similarity[i] for i in idx]
        lineages = [NamedTuple{ranks}(Tuple(hits[i, r] for r in ranks)) for i in idx]
        path = _lca_path(sims, lineages, min_agreement, limits, top_delta, ignore_missing,
                         valid_name, ranks)
        row = Dict{Symbol,Any}(:sequence => s, :n_hits => length(idx))
        if isempty(idx)
            row[:best_similarity] = missing
            row[:best_process_id] = missing
        else
            b = idx[argmax(sims)]
            row[:best_similarity] = hits.similarity[b]
            row[:best_process_id] = hits.process_id[b]
        end
        row[:assigned_rank] = isempty(path) ? missing : String(last(path).rank)
        row[:taxon] = isempty(path) ? missing : last(path).taxon
        for r in ranks
            row[r] = missing
            row[Symbol(r, :_agreement)] = missing
        end
        for p in path
            row[p.rank] = p.taxon
            row[Symbol(p.rank, :_agreement)] = p.agreement
        end
        push!(out, row)
    end
    return out
end

# Rows of each normalized sequence holding a hit, sequences in order of first
# appearance and the number of exact repeats dropped. The same process ID with
# different data for one sequence is an error, since it would be counted twice.
function _hit_groups(hits, ranks)
    groups = Dict{String,Vector{Int}}()
    order = String[]
    first_row = Dict{Tuple{String,String},Int}()
    repeated = 0
    for i in 1:nrow(hits)
        raw = hits.sequence[i]
        ismissing(raw) && throw(ArgumentError("hit table row $i has a missing sequence"))
        s = normalize_sequence(raw)
        haskey(groups, s) || (push!(order, s); groups[s] = Int[])
        sim = hits.similarity[i]
        ismissing(sim) && continue
        (isfinite(sim) && 0 <= sim <= 100) ||
            throw(ArgumentError("hit table row $i has similarity $sim outside [0, 100]"))
        pid = hits.process_id[i]
        if !ismissing(pid)
            key = (s, String(pid))
            j = get(first_row, key, 0)
            if j > 0
                same = isequal(hits.similarity[j], sim) &&
                       all(r -> isequal(hits[j, r], hits[i, r]), ranks)
                same || throw(ArgumentError("process ID $pid appears twice with different data " *
                                            "for the same sequence (rows $j and $i)"))
                repeated += 1
                continue
            end
            first_row[key] = i
        end
        push!(groups[s], i)
    end
    return groups, order, repeated
end
