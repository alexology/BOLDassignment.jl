const _IUPAC = Set("ACGTURYSWKMBDHVN")

"""
    normalize_sequence(seq) -> String

Remove whitespace and alignment gaps (`-`, `.`), convert to upper case and
`U` to `T`. Throws `ArgumentError` for characters outside the IUPAC
nucleotide alphabet or when nothing is left.
"""
function normalize_sequence(seq::AbstractString)
    io = IOBuffer()
    for c in seq
        (isspace(c) || c == '-' || c == '.') && continue
        u = uppercase(c)
        u in _IUPAC || throw(ArgumentError("invalid nucleotide '$c' in sequence"))
        print(io, u == 'U' ? 'T' : u)
    end
    s = String(take!(io))
    isempty(s) && throw(ArgumentError("empty sequence"))
    return s
end

# FASTA payload of one submission; identifiers are the positions in `seqs`.
function _fasta_bytes(seqs::AbstractVector{String}, ids::AbstractVector{String})
    io = IOBuffer()
    writer = FASTA.Writer(io; width = -1)
    for (id, s) in zip(ids, seqs)
        write(writer, FASTA.Record(id, s))
    end
    flush(writer)
    return take!(io)
end
