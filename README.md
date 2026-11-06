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

## Install from GitHub

This package is not registered in the Julia General registry. Install it
from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/alexology/BOLDassignment.jl")
```

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

## Assignment rule

Two directions are available (`direction` keyword of `lca_path`,
`lca_assign` and `bold_assign`). The default, `:top_down`, descends from
phylum to species; `:bottom_up` starts from the species and moves up.

### Top-down (default)

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

### Bottom-up

With `direction = :bottom_up` the ranks are visited from species to phylum,
and the first rank whose pool agrees is the **decision rank**:

- the pool of a rank is defined as above (threshold and `top_delta`);
- only hits with a valid name at that rank vote; with `ignore_missing = true`
  the others leave the pool of that rank, with `ignore_missing = false` they
  stay and count against the agreement. Each rank is judged on its own: a
  record of `Protonemura sp.` does not vote for the species but still votes
  for the genus;
- the most frequent taxon is accepted if it is unique (no tie) and supported
  by at least `min_agreement` percent of the pool; otherwise, or with an empty
  pool, the next more general rank is tried. If no rank agrees, the sequence
  is not assigned.

The ranks above the decision rank are taken from the **supporting hits** (the
pool hits voting for the accepted taxon): at each rank, the most frequent
valid name among them. For these ranks `<rank>_agreement` is the percent of
supporting hits with that name, not a share of the rank's pool. If they
disagree (e.g. the same species name under two genera in BOLD), the most
frequent name is kept and `lineage_conflict` is true; on a tie that rank and
the more general ones stay `missing`. With this direction `lca_assign` and
`bold_assign` add the columns `decision_rank` and `lineage_conflict`; with
the default direction the columns are unchanged.

**Example.** Best hit 99.3 % on *Protonemura intricata* (Plecoptera); 45
hits on Plecoptera between 85.9 and 99.3 %, 55 on Diptera between 85.4 and
86.5 %. Top-down with `top_delta = Inf`, the order pool (≥ 85 %) contains
both orders (45 against 55), so the assignment stops at the class *Insecta*.
Bottom-up, the species pool (≥ 97 %) holds only *P. intricata* hits, so the
species is assigned, and genus, family, order, class and phylum follow from
those hits.

**Which one to use.** Top-down is the strict LCA: a rank is assigned only
when all the hits above its threshold agree on every more general rank, so
distant hits of other groups can stop the assignment early, as in the example.
It suits cautious inventories and incomplete reference libraries, where a
close match may still be misleading. Bottom-up asks, for each rank, whether
the hits close enough for that rank agree; distant hits do not block a
well-supported species. It is closer to how a close match is read by hand and
to threshold-based classifiers. In both, a valid congener above the species
threshold (e.g. 97.5 % against a best hit of 100 %) blocks the species with
`min_agreement = 100`; bottom-up then assigns the genus. Check
`lineage_conflict` before relying on the higher ranks.

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
does not finish within `timeout` is resubmitted at most `max_resubmits` times;
the error then names the last submission. If a finished submission lacks the
result of some sequences, they are reported with a warning and returned
without hits.

With `checkpoint = dir`, the raw result of every sequence is stored in `dir`
(one file per sequence; file name = SHA-256 of database, mode and sequence)
and reused on later calls: an interrupted run resumes without resubmitting
finished sequences, and adding sequences submits only the new ones.
Sequences missing from BOLD results are not stored, so they are submitted
again. Files are written to a unique temporary name and renamed; a
checkpoint that cannot be read back is reported and its sequence submitted
again. Checkpoints written by versions before 2026-10-25 (one file per
chunk) are not reused.

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

## License

MIT (see `LICENSE`). BOLD data retrieved with this package are subject to
the terms of use of the BOLD Systems.
