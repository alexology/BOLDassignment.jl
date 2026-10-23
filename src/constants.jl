"""
    BOLD_DATABASES

Reference libraries accepted by the identification engine (`db` parameter).
Keys are the short names used by this package, values the API codes.
"""
const BOLD_DATABASES = Dict{Symbol,String}(
    :animal_public => "public.tax-derep",
    :animal_all => "all.tax-derep",
    :animal_species => "species",
    :animal_secondary => "all.animal-alt",
    :plants => "public.plants",
    :fungi => "public.fungi",
    :canadian_arthropods => "DS-CANREF22",
    :red_list => "DS-IUCNPUB",
)

"""
    BOLD_MODES

Operating modes documented by the identification engine: minimum identity
`mi`, minimum overlap `mo`, maximum hits per sequence `maxh` and maximum
number of sequences per submission `max_sequences`.
"""
const BOLD_MODES = Dict{Symbol,NamedTuple{(:mi, :mo, :maxh, :max_sequences),Tuple{Float64,Int,Int,Int}}}(
    :rapid => (mi = 0.94, mo = 100, maxh = 25, max_sequences = 1000),
    :genus_species => (mi = 0.90, mo = 100, maxh = 50, max_sequences = 200),
    :exhaustive => (mi = 0.75, mo = 100, maxh = 100, max_sequences = 100),
)

"""
    RANKS

Taxonomic ranks used for the assignment, from the most general to the most
specific.
"""
const RANKS = (:phylum, :class, :order, :family, :genus, :species)

"""
    DEFAULT_THRESHOLDS

Minimum percent similarity of a hit to support each rank. Same values as
BOLDigger3.
"""
const DEFAULT_THRESHOLDS = Dict{Symbol,Float64}(
    :phylum => 50.0, :class => 75.0, :order => 85.0,
    :family => 90.0, :genus => 95.0, :species => 97.0,
)

_database_code(db::Symbol) = haskey(BOLD_DATABASES, db) ? BOLD_DATABASES[db] :
    throw(ArgumentError("unknown database :$db; valid: $(join(sort!(collect(keys(BOLD_DATABASES))), ", "))"))
_database_code(db::AbstractString) = String(db)

_mode(mode::Symbol) = haskey(BOLD_MODES, mode) ? BOLD_MODES[mode] :
    throw(ArgumentError("unknown mode :$mode; valid: $(join(sort!(collect(keys(BOLD_MODES))), ", "))"))
