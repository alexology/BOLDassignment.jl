# BOLDassignment.jl

Taxonomic assignment of amplicon sequence variants (ASVs) against the
[BOLD v5 identification engine](https://id.boldsystems.org), with a lowest
common ancestor (LCA) that accepts a minimum agreement percentage.

The package uses the documented Barcode ID API
(`https://id.boldsystems.org/openapi.json`):

1. `POST /submission` with a FASTA file and the database/mode parameters;
2. `GET /submission/status/{sub_id}` until the submission is complete;
3. `GET /submission/results/{sub_id}` for the hits (JSON Lines).

Record metadata (BCDM fields) come from the portal API
(`https://portal.boldsystems.org/openapi.json`).

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

`bold_assign` keeps every row of the input, in order, and appends `n_hits`,
`best_similarity`, `best_process_id`, `assigned_rank`, `taxon`, one column per
rank (`phylum` … `species`) and one `<rank>_agreement` column per rank. Use
`prefix = "bold_"` when the input already has columns with these names.

## Assignment rule

For each sequence and each rank, from phylum to species:

- the **pool** contains the hits with similarity ≥ the rank threshold
  (`thresholds`, default as BOLDigger3: phylum 50, class 75, order 85,
  family 90, genus 95, species 97 %) and ≥ best similarity − `top_delta`
  (default `Inf`, i.e. thresholds only);
- the most frequent taxon among the pool hits that agree with the ranks
  already assigned is accepted when it is supported by at least
  `min_agreement` percent of the pool (default **100**, strict LCA);
- otherwise, or on a tie for the most frequent taxon, the assignment stops.

`top_delta = Inf` is a deliberate default: every hit above the rank
threshold takes part, so a congener at 97.5 % keeps the species open even
when the best hit is at 100 %. This is the conservative choice; a finite
value (e.g. `top_delta = 1`) lets the closest hits decide, closer to
top-hit methods such as BOLDigger3. Neither setting handles misidentified
records close to the best hit: use `min_agreement` below 100 for those.

Only hits **identified at a rank** vote at that rank. A hit whose name is
missing or not valid (`is_valid_name`) leaves the pool at that rank and at the
deeper ones. At species level this removes open nomenclature and records
identified only to genus or above: `Onthophagus sp.`, `Onthophagus sp. 1`,
`Onthophagus cf. taurus`, `aff.`, `nr.`, `complex`, interim codes such as
`Diptera sp. BIOUG1234-A01`. Such a hit above 97 % therefore neither blocks
nor supports the species, but still supports the genus. With
`ignore_missing = false` these hits stay in the pool and count against the
agreement; `valid_name = (rank, name) -> ...` replaces the name check.
`lca_path` applies the rule to a single set of hits.

## Databases and modes

| `db` | API code |
|---|---|
| `:animal_public` (default) | `public.tax-derep` |
| `:animal_all` | `all.tax-derep` |
| `:animal_species` | `species` |
| `:animal_secondary` | `all.animal-alt` |
| `:plants` | `public.plants` |
| `:fungi` | `public.fungi` |
| `:canadian_arthropods` | `DS-CANREF22` |
| `:red_list` | `DS-IUCNPUB` |

| `mode` | `mi` | `maxh` | max sequences per submission |
|---|---:|---:|---:|
| `:rapid` (default) | 0.94 | 25 | 1000 |
| `:genus_species` | 0.90 | 50 | 200 |
| `:exhaustive` | 0.75 | 100 | 100 |

`mo` (minimum overlap) is 100 in all modes. A sequence is assigned only from
the hits returned by the chosen mode: in `:rapid` no hit is below 94 %.

## Rate limits

The API documents at most **3 submissions/minute** and **40 status/result
requests/minute**; exceeding them may block the client. `IDEngine()` waits
21 s between submissions and 1.6 s between other requests, keeps at most 5
submissions open and retries failures (HTTP 408, 429, 5xx, network errors)
with increasing waits. All settings are keywords of `IDEngine`.

The spacing is shared by every engine with the default `throttle`, so
consecutive or concurrent calls respect the limits together. A submission that
does not finish within `timeout` is resubmitted at most `max_resubmits` times.

With `checkpoint = dir`, the raw response of every chunk is stored in `dir`
(file name = SHA-256 of database, mode and sequences) and reused on later
calls, so an interrupted run resumes without resubmitting finished chunks.
Sequences are sorted before chunking, so the input order does not matter.
Files are written to a unique temporary name and renamed; a checkpoint that
cannot be read back is reported and its chunk submitted again.

## Input checks

- `bold_assign` checks every argument (LCA options, thresholds, ranks,
  column names, sequences) **before** contacting BOLD.
- `thresholds` may list only the ranks to change, e.g.
  `Dict(:species => 98)`; the other ranks keep the defaults. Unknown ranks
  and values outside [0, 100] are errors.
- Repeated hit rows (same sequence and process ID, same data, e.g. after
  `vcat` of two tables) are counted once with a warning; the same process
  ID with different data is an error.
- Hit sequences are normalized like the input, so case and gaps do not matter.

## Tests

```
julia --project=. test/runtests.jl
```

The suite uses a simulated engine and makes no network request. Set
`BOLDASSIGNMENT_LIVE=1` to add a test against the real BOLD services.
