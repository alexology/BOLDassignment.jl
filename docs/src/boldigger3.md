# Top hits as in BOLDigger3

## BOLDigger3 top hits

`boldigger3_select` is a translation of the top-hit selection of
[BOLDigger3](https://github.com/DominikBuchner/BOLDigger3) 3.0.2
(`select_top_hit.py`, MIT licence, see `NOTICE`), with BOLDigger3's output
columns (`records`, `records_ratio`, `selected_level`, `BIN`, `flags`).
The rule differs from the LCA above: the best similarity picks the starting
level (97/95/90/85/75 %, plus 50 % for phyla), the most numerous taxon among
the hits at or above that threshold wins and the flags report its weak points.
Names with digits or punctuation other than `-` (e.g. `Baetis sp. 1`) are
treated as missing, as in BOLDigger3.

```julia
hits = bold_identify(asvs.ASV; db = :animal_public, mode = :rapid)
ids = Dict(zip(normalize_sequence.(asvs.ASV), asvs.id))
top = boldigger3_select(boldigger3_hits(hits; ids, metadata = meta))
```

BOLDigger3 takes the record status (public or private) and the
identification method (flag 1) from its local database (`bold_public.*.ddb`):
pass them as `metadata` (columns `processid`, `identification_method`, one
row per record and marker, as in that database) to obtain the same flags.
`validation/boldigger3_reference.py metadata` extracts them.

**Ties.** When several taxa share the highest record count, BOLDigger3 takes
the first one after a pandas sort that is not stable (numpy's SIMD sort, so
its choice can depend on the numpy version and the CPU); `boldigger3_select`
takes the first one in hit order and reports the number of tied taxa in
`tied_groups`. BOLDigger3 also receives hits of equal similarity in an order
that its DuckDB query does not define, which can change the order of the
BINs in `BIN` and the winner of a tie.

**Validation** (`validation/boldigger3_*.{py,jl}`, BOLDigger3 3.0.2, pandas
2.2.3, numpy 2.3.2, DuckDB 1.3.2, 2026-11-14):

- *Same input* (the hit table BOLDigger3's `find_top_hit` receives): 1,303 of
  1,303 real COI ASVs (Element AVITI, BF2–BR2) identical in every column,
  including 11 ties; 2,000 random hit tables (thresholds, missing and invalid
  names, no-matches, private records, methods, BINs) identical in all 1,559
  without ties, and in 436 of the 441 ties. `test/data/boldigger3` holds 600
  of these cases for the test suite.
- *End to end* (`bold_identify` + metadata from the same `.ddb` +
  `boldigger3_hits` + `boldigger3_select`, against a separate BOLDigger3 run on
  the same ASVs): same BOLD records and similarities for all 1,303 ASVs;
  1,289 of the 1,292 ASVs without ties identical in every column (the other 3
  differ only in the order of the BINs), 8 of the 11 ties with the same taxon.
