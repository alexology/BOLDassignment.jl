const _LCA_KEYWORDS = (:min_agreement, :thresholds, :top_delta, :ignore_missing, :valid_name, :ranks)

"""
    bold_assign(df::AbstractDataFrame; sequence_col = :ASV, hits = nothing, prefix = "",
                db = :animal_public, mode = :rapid, engine = IDEngine(),
                chunk_size = nothing, checkpoint = nothing, lca_kwargs...) -> DataFrame

Assign a taxonomy to the sequences in column `sequence_col` of `df`.

Unique sequences are identified with [`bold_identify`](@ref) (or taken from a
previously computed `hits` table, without network access) and assigned with
[`lca_assign`](@ref); `lca_kwargs` (`min_agreement`, `thresholds`,
`top_delta`, `ignore_missing`, `valid_name`, `ranks`) are passed to it. Every row of `df`
is kept, in order, with the assignment columns appended, their names preceded
by `prefix`. Rows with a `missing` sequence get `missing` assignments.

# Example
```julia
asvs = DataFrame(ASV = ["ACTTTATATTTTATTTTTGGAATTTGAGCAGGAATAGTAGG..."])
hits = bold_identify(asvs.ASV; checkpoint = "bold_cache")
tax = bold_assign(asvs; hits, min_agreement = 90)
```
"""
function bold_assign(df::AbstractDataFrame; sequence_col = :ASV, hits = nothing,
                     prefix::AbstractString = "", db = :animal_public, mode::Symbol = :rapid,
                     engine::IDEngine = IDEngine(), chunk_size = nothing, checkpoint = nothing,
                     lca_kwargs...)
    for k in keys(lca_kwargs)
        k in _LCA_KEYWORDS || throw(ArgumentError("unknown keyword $k"))
    end
    # Every argument is checked before the (possibly long) identification.
    ranks = _check_ranks(get(lca_kwargs, :ranks, RANKS))
    _check_thresholds(get(lca_kwargs, :thresholds, DEFAULT_THRESHOLDS), ranks)
    _check_lca_options(get(lca_kwargs, :min_agreement, 100), get(lca_kwargs, :top_delta, Inf))
    col = Symbol(sequence_col)
    hasproperty(df, col) || throw(ArgumentError("column :$col not found"))
    new_columns = [Symbol(prefix, n) for n in (_ASSIGNMENT_COLUMNS[2:end]..., ranks...,
                                                (Symbol(r, :_agreement) for r in ranks)...)]
    for n in new_columns
        hasproperty(df, n) &&
            throw(ArgumentError("column :$n already exists; use the `prefix` keyword"))
    end
    normalized = map(enumerate(df[!, col])) do (i, s)
        ismissing(s) && return missing
        try
            normalize_sequence(s)
        catch err
            err isa ArgumentError || rethrow()
            throw(ArgumentError("row $i of :$col: $(err.msg)"))
        end
    end
    if hits === nothing
        hits = bold_identify(collect(skipmissing(normalized)); db, mode, engine, chunk_size,
                             checkpoint)
    end
    table = lca_assign(hits; lca_kwargs...)
    position = Dict(s => i for (i, s) in enumerate(table.sequence))
    for (i, s) in enumerate(normalized)
        ismissing(s) || haskey(position, s) ||
            throw(ArgumentError("hits table has no entry for row $i of :$col"))
    end
    out = copy(df)
    for name in propertynames(table)
        name === :sequence && continue
        newname = Symbol(prefix, name)
        src = table[!, name]
        out[!, newname] = [ismissing(s) ? missing : src[position[s]] for s in normalized]
    end
    return out
end
