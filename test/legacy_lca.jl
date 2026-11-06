# Top-down LCA of BOLDassignment at commit 63e5d70 (version 0.1.2), copied
# verbatim by a script except for the function names and the module prefix of
# internal helpers. Used only to check that the default direction did not
# change; do not edit.

function legacy_lca_path(similarity, lineages, min_agreement, limits, top_delta, ignore_missing,
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

function legacy_lca_assign(hits::AbstractDataFrame; min_agreement::Real = 100,
                    thresholds = BOLDassignment.DEFAULT_THRESHOLDS, top_delta::Real = Inf,
                    ignore_missing::Bool = true, valid_name = BOLDassignment.is_valid_name,
                    ranks = BOLDassignment.RANKS)
    ranks = BOLDassignment._check_ranks(ranks)
    limits = BOLDassignment._check_thresholds(thresholds, ranks)
    BOLDassignment._check_lca_options(min_agreement, top_delta)
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
    groups, order, repeated = BOLDassignment._hit_groups(hits, ranks)
    repeated > 0 && @warn "Ignored $repeated repeated hit row(s) (same sequence, process ID and data)"
    for s in order
        idx = groups[s]
        sims = Float64[hits.similarity[i] for i in idx]
        lineages = [NamedTuple{ranks}(Tuple(hits[i, r] for r in ranks)) for i in idx]
        path = legacy_lca_path(sims, lineages, min_agreement, limits, top_delta, ignore_missing,
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
