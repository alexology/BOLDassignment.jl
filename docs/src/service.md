# The BOLD service

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
