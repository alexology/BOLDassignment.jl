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

const LCA_DIRECTIONS = (:top_down, :bottom_up)

function _check_direction(direction)
    direction isa Symbol && direction in LCA_DIRECTIONS ||
        throw(ArgumentError("direction must be :top_down or :bottom_up, got $(repr(direction))"))
    return direction
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
             ranks = RANKS, direction = :top_down)

Lowest common ancestor of the hits of one sequence with a minimum agreement.
`direction = :top_down` (default) descends from the most general rank, as
described first; `direction = :bottom_up` starts from the most specific rank
and moves up (see below).

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
one per assigned rank from the most general to the most specific (empty when
nothing can be assigned).

# Bottom-up direction

With `direction = :bottom_up` the ranks are visited from the most specific to
the most general, and the **decision rank** is the first one whose pool
agrees. The pool of a rank is defined as above. Only hits with a valid name at
that rank vote; with `ignore_missing = true` the others leave the pool of that
rank, with `ignore_missing = false` they stay in it and count against the
agreement. Unlike the top-down direction, every rank is judged on its own: a
hit not identified at the species still votes at the genus, and the votes at
a rank do not depend on the other ranks. The most frequent taxon is accepted
when it is unique and supported by at least `min_agreement` percent of the
pool; otherwise (or with an empty pool) the next more general rank is tried.
When no rank agrees, nothing is assigned.

The ranks above the decision rank come from the **supporting hits**, the pool
hits voting for the accepted taxon: at each higher rank the most frequent
valid name among them is taken, even when they are not unanimous. In those
entries `n_pool` is the number of supporting hits, `n_support` the number
carrying the name and `agreement = 100 * n_support / n_pool`, so an agreement
below 100 signals inconsistent BOLD classifications (e.g. one species name
under two genera, or supporting hits without a name at that rank). On a tie,
or when no supporting hit has a valid name, that rank and all the more general
ones are left out of the path. The decision rank is the last entry and the
only one to which `min_agreement` applies. Use `top_delta = Inf` (default)
with this direction, so that every hit above the threshold of a rank takes
part.
"""
function lca_path(similarity::AbstractVector{<:Real}, lineages::AbstractVector;
                  min_agreement::Real = 100, thresholds = DEFAULT_THRESHOLDS,
                  top_delta::Real = Inf, ignore_missing::Bool = true,
                  valid_name = is_valid_name, ranks = RANKS, direction = :top_down)
    _check_direction(direction)
    ranks = _check_ranks(ranks)
    limits = _check_thresholds(thresholds, ranks)
    _check_lca_options(min_agreement, top_delta)
    length(similarity) == length(lineages) ||
        throw(ArgumentError("similarity and lineages must have the same length"))
    for s in similarity
        (isfinite(s) && 0 <= s <= 100) ||
            throw(ArgumentError("similarity must be finite and within [0, 100], got $s"))
    end
    direction === :bottom_up &&
        return first(_lca_path_bottom_up(similarity, lineages, min_agreement, limits, top_delta,
                                         ignore_missing, valid_name, ranks))
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

const _PathEntry = NamedTuple{(:rank, :taxon, :agreement, :n_support, :n_pool),
                              Tuple{Symbol,String,Float64,Int,Int}}

# Bottom-up rule of `lca_path`. Returns the path and whether the supporting
# hits disagree (or lack a valid name) at a rank above the decision rank.
function _lca_path_bottom_up(similarity, lineages, min_agreement, limits, top_delta,
                             ignore_missing, valid_name, ranks)
    path = _PathEntry[]
    isempty(similarity) && return path, false
    floor_sim = maximum(similarity) - top_delta
    counts = Dict{String,Int}()
    voters = Int[]
    for j in reverse(eachindex(ranks))
        r = ranks[j]
        t = limits[r]
        empty!(counts)
        empty!(voters)
        n_pool = 0
        for i in eachindex(similarity)
            s = similarity[i]
            (s >= t && s >= floor_sim) || continue
            name = lineages[i][r]
            if valid_name(r, name)
                n_pool += 1
                push!(voters, i)
                counts[name] = get(counts, name, 0) + 1
            elseif !ignore_missing
                n_pool += 1
            end
        end
        (n_pool == 0 || isempty(counts)) && continue
        best = maximum(values(counts))
        count(==(best), values(counts)) == 1 || continue
        taxon = first(k for (k, v) in counts if v == best)
        agreement = 100 * best / n_pool
        # Same relative tolerance as the top-down rule.
        agreement >= min_agreement * (1 - 1e-12) || continue
        supporting = [i for i in voters if lineages[i][r] == taxon]
        above, conflict = _ranks_above(lineages, supporting, ranks, j, valid_name)
        append!(path, above)
        push!(path, (rank = r, taxon = taxon, agreement = agreement, n_support = best,
                     n_pool = n_pool))
        return path, conflict
    end
    return path, false
end

# Names of the ranks above `ranks[j]` taken from the supporting hits, most
# general first. Stops at the first tie or rank without any valid name.
function _ranks_above(lineages, supporting, ranks, j, valid_name)
    entries = _PathEntry[]
    conflict = false
    n = length(supporting)
    counts = Dict{String,Int}()
    for q in j-1:-1:1
        r = ranks[q]
        empty!(counts)
        for i in supporting
            name = lineages[i][r]
            valid_name(r, name) && (counts[name] = get(counts, name, 0) + 1)
        end
        if isempty(counts)
            conflict = true
            break
        end
        best = maximum(values(counts))
        if count(==(best), values(counts)) != 1
            conflict = true
            break
        end
        best == n || (conflict = true)
        taxon = first(k for (k, v) in counts if v == best)
        pushfirst!(entries, (rank = r, taxon = taxon, agreement = 100 * best / n,
                             n_support = best, n_pool = n))
    end
    return entries, conflict
end

"""
    lca_assign(hits::DataFrame; min_agreement = 100, thresholds = DEFAULT_THRESHOLDS,
               top_delta = Inf, ignore_missing = true, valid_name = is_valid_name,
               ranks = RANKS, direction = :top_down) -> DataFrame

Apply [`lca_path`](@ref) to every sequence of a hit table returned by
[`bold_identify`](@ref).

One row per sequence, in order of first appearance, with columns `sequence`,
`n_hits`, `best_similarity`, `best_process_id`, `assigned_rank`, `taxon`
(deepest assigned taxon), one column per rank with the assigned name and one
`<rank>_agreement` column with the percent of pool hits supporting it.
Unassigned ranks are `missing`; a sequence without hits has `n_hits = 0`.

With `direction = :bottom_up` (see [`lca_path`](@ref)) two columns follow:
`decision_rank`, the rank decided by its own pool (equal to `assigned_rank`,
`missing` when nothing is assigned), and `lineage_conflict`, true when the
supporting hits disagree or lack a valid name at a rank above it. For the
ranks above the decision rank, `<rank>_agreement` is the percent of
supporting hits carrying that name, not a share of the rank's pool. With the
default `:top_down` the columns are exactly as in earlier versions.
"""
function lca_assign(hits::AbstractDataFrame; min_agreement::Real = 100,
                    thresholds = DEFAULT_THRESHOLDS, top_delta::Real = Inf,
                    ignore_missing::Bool = true, valid_name = is_valid_name,
                    ranks = RANKS, direction = :top_down)
    _check_direction(direction)
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
    bottom_up = direction === :bottom_up
    if bottom_up
        out[!, :decision_rank] = Union{Missing,String}[]
        out[!, :lineage_conflict] = Bool[]
    end
    groups, order, repeated = _hit_groups(hits, ranks)
    repeated > 0 && @warn "Ignored $repeated repeated hit row(s) (same sequence, process ID and data)"
    for s in order
        idx = groups[s]
        sims = Float64[hits.similarity[i] for i in idx]
        lineages = [NamedTuple{ranks}(Tuple(hits[i, r] for r in ranks)) for i in idx]
        path, conflict = if bottom_up
            _lca_path_bottom_up(sims, lineages, min_agreement, limits, top_delta,
                                ignore_missing, valid_name, ranks)
        else
            _lca_path(sims, lineages, min_agreement, limits, top_delta, ignore_missing,
                      valid_name, ranks), false
        end
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
        if bottom_up
            row[:decision_rank] = row[:assigned_rank]
            row[:lineage_conflict] = conflict
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
