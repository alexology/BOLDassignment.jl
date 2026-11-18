# The LCA assignment

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
