# BOLDassignment.jl

Taxonomic assignment of amplicon sequence variants (ASVs) against the
[BOLD v5 identification engine](https://id.boldsystems.org), with a lowest
common ancestor (LCA) that accepts a minimum agreement percentage.

This package is part of an experimental effort to use AI to port
bioinformatics and quantitative ecology packages to Julia, with the goal of
creating a unified environment for high performance scientific computing.

The package uses the documented Barcode ID API
(`https://id.boldsystems.org/openapi.json`):

1. `POST /submission` with a FASTA file and the database/mode parameters;
2. `GET /submission/status/{sub_id}` until the submission is complete;
3. `GET /submission/results/{sub_id}` for the hits (JSON Lines).

Record metadata (BCDM fields) come from the portal API
(`https://portal.boldsystems.org/openapi.json`).

## Installation

This package is not registered in the Julia General registry. Install it
from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/alexology/BOLDassignment.jl")
```

The development versions published on GitHub were numbered 0.1.1–0.1.3;
the numbering restarts at 0.1.0 for the first registration (the General
registry accepts only 0.0.1, 0.1.0 or 1.0.0 as a first version). The
functionality is unchanged.

## Usage

```julia
using BOLDassignment, DataFrames

asvs = DataFrame(id = ["asv1", "asv2"], ASV = [seq1, seq2])

# One call: identification + strict LCA (100% agreement).
tax = bold_assign(asvs)

# Two steps: keep the hits and try different LCA settings offline.
hits = bold_identify(asvs.ASV; db = :animal_public, mode = :rapid, checkpoint = "bold_cache")
strict   = bold_assign(asvs; hits)
majority = bold_assign(asvs; hits, min_agreement = 90)

# Collection data of the best hits.
meta = bold_records(skipmissing(majority.best_process_id) |> collect;
                    fields = [:processid, :bin_uri, Symbol("country/ocean"), :coord])
```

The hit table of `bold_identify` has `sequence`, `process_id`, `bin_uri`,
`record_status`, `similarity` (`100 - pdist`), `marker` (e.g. `COI-5P`),
`taxid`, `align_length` (`length` of the alignment), `query_coverage`
(`qcovhsp`, or computed from `qstart`/`qend`), `result_received` and one
column per rank. `result_received = false` marks a sequence whose result
BOLD did not return in a finished submission: its hits are unknown (not
"no match") and it is submitted again by a later call with a checkpoint.

The portal returns one record per specimen **and marker** (`record_id`
`<processid>.<marker>`, each with its own sequence `nuc`): use
`bold_records(ids; marker = "COI-5P")` to keep only the records of one marker.

`bold_assign` keeps every row of the input, in order, and appends `n_hits`,
`best_similarity`, `best_process_id`, `assigned_rank`, `taxon`, one column per
rank (`phylum` … `species`) and one `<rank>_agreement` column per rank. Use
`prefix = "bold_"` when the input already has columns with these names.

## Tests

```
julia --project=. test/runtests.jl
```

The suite uses a simulated engine and makes no network request. Set
`BOLDASSIGNMENT_LIVE=1` to add a test against the real BOLD services.

## License

MIT (see `LICENSE`). `src/boldigger3.jl` translates code of BOLDigger3
(MIT, Copyright (c) 2024 DominikBuchner; see `NOTICE`). BOLD data retrieved with this package are subject to
the terms of use of the BOLD Systems.
