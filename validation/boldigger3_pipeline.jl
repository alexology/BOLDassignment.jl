# End-to-end comparison with BOLDigger3: bold_identify hits (saved as TSV)
# + metadata from BOLDigger3's local database -> boldigger3_hits ->
# boldigger3_select, against BOLDigger3's own result for the same sequences.
#   julia --project=. validation/boldigger3_pipeline.jl <bold_hits.tsv> <asvs.fasta> \
#         <metadata.tsv> <reference dir (boldigger3_reference.py real)>
# The BOLD hits of the two tools come from separate requests: the script
# first reports, per sequence, whether they are the same records with the
# same similarities, then the agreement of the assignments.
using BOLDassignment, DataFrames

include(joinpath(@__DIR__, "boldigger3_compare.jl"))   # read_tsv, num

function main(hits_path, fasta, meta_path, refdir)
    hits = read_tsv(hits_path)
    hits.similarity = num.(hits.similarity)
    ids = Dict{String,String}(); name = ""
    for l in eachline(fasta)
        startswith(l, '>') ? (name = split(l[2:end])[1]) : (ids[uppercase(strip(l))] = name)
    end
    meta = read_tsv(meta_path)
    ref_in = read_tsv(joinpath(refdir, "input_hits.tsv"))
    ref_in.pct_identity = num.(ref_in.pct_identity)
    exp = read_tsv(joinpath(refdir, "expected.tsv"))

    # 1. the BOLD records of each sequence
    records(df, idcol, pcol) = Dict(k => Set(zip(g.process_id, g[!, pcol])) for (k, g) in
                                    pairs(groupby(dropmissing(df, :process_id), idcol)) for k in (k[1],))
    ours = records(transform(hits, :sequence => ByRow(s -> ids[s]) => :id), :id, :similarity)
    theirs = records(ref_in, :id, :pct_identity)
    same_hits = Set(i for i in exp.id if get(ours, i, Set()) == get(theirs, i, Set()))
    println("BOLD records and similarities identical: ", length(same_hits), "/", nrow(exp), " sequences")

    # 2. assignments
    got = boldigger3_select(boldigger3_hits(hits; ids, metadata = meta))
    got = got[indexin(exp.id, got.id), :]
    tied = parse.(Int, exp.top_ties) .> 1
    cols = (RANKS..., :selected_level, :records, :records_ratio, :BIN, :flags, :status, :pct_identity)
    eq(c) = c in (:pct_identity, :records_ratio) ? isequal.(got[!, c], num.(exp[!, c])) :
            c == :records ? isequal.(got[!, c], parse.(Int, exp[!, c])) : isequal.(got[!, c], exp[!, c])
    sameset = [i in same_hits for i in exp.id]
    rows = trues(nrow(exp))
    for c in cols
        s = eq(c); rows .&= s
        println(rpad(string(c), 16), "identical: ", count(s[sameset .& .!tied]), "/", count(sameset .& .!tied),
                " (same records, no tie), ", count(s[sameset .& tied]), "/", count(sameset .& tied),
                " (same records, tie), ", count(s[.!sameset]), "/", count(.!sameset), " (different records)")
    end
    println("whole rows identical: ", count(rows[sameset .& .!tied]), "/", count(sameset .& .!tied),
            " (same records, no tie)")
    for i in findall(.!rows .& sameset .& .!tied)[1:min(end, 10)]
        println("  ", exp.id[i], ": ", join(("$c Julia $(repr(got[i, c])) BOLDigger3 $(repr(exp[i, c]))"
                                              for c in cols if !eq(c)[i]), "; "))
    end
    return (; got, exp, tied, sameset, rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS...)
end
