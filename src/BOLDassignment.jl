"""
    BOLDassignment

Taxonomic assignment of amplicon sequence variants (ASVs) against the BOLD v5
identification engine (<https://id.boldsystems.org>).

Sequences are submitted through the documented Barcode ID API, the hits are
returned as a long table ([`bold_identify`](@ref)) and a lowest common
ancestor with a minimum agreement percentage turns them into one assignment
per sequence ([`lca_assign`](@ref)). [`bold_assign`](@ref) runs both steps on
a `DataFrame` with an `ASV` column. [`bold_records`](@ref) retrieves record
metadata from the BOLD portal API.
"""
module BOLDassignment

using DataFrames
using FASTX
using HTTP
using JSON3
using SHA

export IDEngine, BOLD_DATABASES, BOLD_MODES, DEFAULT_THRESHOLDS, RANKS,
       bold_identify, lca_assign, lca_path, is_valid_name, bold_assign, bold_records,
       normalize_sequence

include("constants.jl")
include("sequences.jl")
include("client.jl")
include("identify.jl")
include("lca.jl")
include("assign.jl")
include("portal.jl")

end
