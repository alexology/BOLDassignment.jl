# Compare boldigger3_select with BOLDigger3 on the output of
# boldigger3_reference.py (input_hits.tsv, expected.tsv).
#   julia --project=. validation/boldigger3_compare.jl <dir>
# Prints the agreement column by column, separately for sequences where
# several taxa share the highest record count (BOLDigger3's choice then
# depends on numpy's sort), and lists the differences.
using BOLDassignment, DataFrames

function read_tsv(path)
    lines = readlines(path)
    header = split(lines[1], '\t')
    cols = [Union{Missing,String}[] for _ in header]
    for l in lines[2:end]
        f = split(l, '\t'; keepempty = true)
        for (j, x) in enumerate(f)
            push!(cols[j], isempty(x) ? missing : x == "\"\"" ? "" : String(x))
        end
    end
    DataFrame([Symbol(h) => c for (h, c) in zip(header, cols)])
end

num(x) = x === missing ? missing : parse(Float64, x)

function compare(dir; verbose = true)
    hits = read_tsv(joinpath(dir, "input_hits.tsv"))
    hits.pct_identity = num.(hits.pct_identity)
    exp = read_tsv(joinpath(dir, "expected.tsv"))
    got = boldigger3_select(hits)
    nrow(got) == nrow(exp) || error("rows: Julia $(nrow(got)), BOLDigger3 $(nrow(exp))")
    all(got.id .== exp.id) || error("sequence order differs")
    tied = parse.(Int, exp.top_ties) .> 1
    same(c) = c in ("pct_identity", "records_ratio") ? isequal.(got[!, c], num.(exp[!, c])) :
              c == "records" ? isequal.(got[!, c], parse.(Int, coalesce.(exp[!, c], "-1"))) :
              isequal.(got[!, c], exp[!, c])
    cols = ["phylum", "class", "order", "family", "genus", "species", "pct_identity", "status",
            "records", "records_ratio", "selected_level", "BIN", "flags"]
    rowsame = trues(nrow(exp))
    for c in cols
        s = same(c); rowsame .&= s
        verbose && println(rpad(c, 16), "identical: ", count(s[.!tied]), "/", count(.!tied),
                           " without ties, ", count(s[tied]), "/", count(tied), " with ties")
    end
    println("whole rows identical: ", count(rowsame[.!tied]), "/", count(.!tied), " without ties, ",
            count(rowsame[tied]), "/", count(tied), " with ties; tied_groups = top_ties: ",
            all(got.tied_groups .== parse.(Int, exp.top_ties)))
    for i in findall(.!rowsame .& .!tied)[1:min(end, 10)]
        println("  differs (no tie): ", exp.id[i])
        for c in cols
            same(c)[i] || println("    ", c, ": Julia ", repr(got[i, c]), ", BOLDigger3 ", repr(exp[i, c]))
        end
    end
    return (; got, exp, tied, rowsame)
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare(ARGS[1])
end
