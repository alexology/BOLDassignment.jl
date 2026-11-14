"""Reference outputs of BOLDigger3's top-hit selection, for the Julia port.

For every sequence this script builds exactly the table that BOLDigger3's
`find_top_hit` receives (same SQL query, then `clean_dataframe`), stores it,
and stores `find_top_hit`'s answer. It also records, for every sequence, how
many groups share the highest record count at the selected level: on such
ties BOLDigger3's choice depends on numpy's (non-stable) argsort.

  python boldigger3_reference.py real <id engine .duckdb> <fasta> <outdir>
      sequences identified by `boldigger3 identify` (final_results table)
  python boldigger3_reference.py metadata <bold_public .ddb> <hits.tsv> <out.tsv>
      processid and identification_method of the process IDs in a hit table
      (column process_id), from BOLDigger3's local database: the `metadata`
      argument of boldigger3_hits for an end-to-end comparison
  python boldigger3_reference.py synthetic <outdir> [n] [seed]
      random hit tables covering thresholds, missing and invalid names,
      no-matches, private records, identification methods and BINs

Outputs (tab-separated, empty field = missing): input_hits.tsv (one row per
hit, in BOLDigger3's order), expected.tsv (one row per sequence) and
versions.txt; `""` marks an empty string. Requires the boldigger3 package (tested with 3.0.2).
"""
import random
import sys
from importlib.metadata import version
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from boldigger3.select_top_hit import clean_dataframe, find_top_hit, get_threshold, move_threshold_up

THRESHOLDS = [97, 95, 90, 85, 75, 50]  # BOLDigger3 defaults (+ virtual 50 for phyla)
LEVELS = ["phylum", "class", "order", "family", "genus", "species"]
INPUT_COLUMNS = ["id", *LEVELS, "pct_identity", "status", "bin_uri", "identification_method", "process_id"]
OUTPUT_COLUMNS = ["id", *LEVELS, "pct_identity", "status", "records", "records_ratio",
                  "selected_level", "BIN", "flags"]
METADATA = ["processid", "sex", "life_stage", "inst", "country/ocean", "identified_by",
            "identification_method", "coord", "nuc", "marker_code"]


def top_count_ties(hits):
    """Groups sharing the highest count at the level find_top_hit selects (its loop, unchanged)."""
    threshold, level = get_threshold(hits, THRESHOLDS)
    if level == "no-match":
        return 0
    while True:
        above = hits.loc[hits["pct_identity"] >= threshold]
        levels = LEVELS[: LEVELS.index(level) + 1]
        groups = (above[levels].groupby(by=levels, sort=False, dropna=False)
                  .size().reset_index(name="count").dropna(subset=level, axis=0))
        if len(groups.index) == 0:
            if level == "phylum":
                return 0
            threshold, level = move_threshold_up(threshold, THRESHOLDS)
            continue
        return int((groups["count"] == groups["count"].max()).sum())


def cell(x):
    if x is None or x is pd.NA or (isinstance(x, float) and np.isnan(x)):
        return ""
    if isinstance(x, (float, np.floating)):
        return repr(float(x))
    return '""' if str(x) == "" else str(x)  # an empty string, unlike a missing value


def write(rows, columns, path):
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write("\t".join(cell(r.get(c)) for c in columns) + "\n")


def run(tables, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    inputs, expected = [], []
    for hits in tables:
        top = find_top_hit(hits.copy(), THRESHOLDS)
        row = {c: (top[c].iloc[0] if c in top.columns else None) for c in OUTPUT_COLUMNS}
        row["top_ties"] = top_count_ties(hits)
        expected.append(row)
        for k, (_, h) in enumerate(hits.iterrows(), start=1):
            r = {c: h[c] for c in INPUT_COLUMNS}
            r["hit"] = k
            inputs.append(r)
    write(inputs, ["id", "hit", *INPUT_COLUMNS[1:]], outdir / "input_hits.tsv")
    write(expected, [*OUTPUT_COLUMNS, "top_ties"], outdir / "expected.tsv")
    with open(outdir / "versions.txt", "w") as f:
        for p in ("boldigger3", "pandas", "numpy", "duckdb"):
            f.write(f"{p} {version(p)}\n")
        f.write(f"thresholds {THRESHOLDS}\n")
    print(f"{len(tables)} sequences, {len(inputs)} hits, "
          f"{sum(r['top_ties'] > 1 for r in expected)} with tied top counts -> {outdir}")


def real_tables(db, fasta):
    ids = [line[1:].split()[0] for line in open(fasta) if line.startswith(">")]
    with duckdb.connect(db, read_only=True) as con:
        for i in ids:
            # BOLDigger3's own query (gather_top_hits), then its cleaning step
            q = f"SELECT * FROM final_results WHERE id='{i}' ORDER BY fasta_order ASC, pct_identity DESC"
            yield clean_dataframe(con.execute(q).df())


# ----------------------------------------------------------------- synthetic --

TREE = [  # phylum, class, order, family, genus, species
    ("Arthropoda", "Insecta", "Ephemeroptera", "Baetidae", "Baetis", "Baetis rhodani"),
    ("Arthropoda", "Insecta", "Ephemeroptera", "Baetidae", "Baetis", "Baetis vardarensis"),
    ("Arthropoda", "Insecta", "Ephemeroptera", "Baetidae", "Cloeon", "Cloeon dipterum"),
    ("Arthropoda", "Insecta", "Diptera", "Chironomidae", "Chironomus", "Chironomus riparius"),
    ("Arthropoda", "Insecta", "Diptera", "Chironomidae", "Polypedilum", "Polypedilum nubeculosum"),
    ("Arthropoda", "Insecta", "Diptera", "Simuliidae", "Simulium", "Simulium ornatum"),
    ("Arthropoda", "Malacostraca", "Amphipoda", "Gammaridae", "Gammarus", "Gammarus fossarum"),
    ("Arthropoda", "Malacostraca", "Amphipoda", "Gammaridae", "Gammarus", "Gammarus pulex"),
    ("Annelida", "Clitellata", "Haplotaxida", "Naididae", "Tubifex", "Tubifex tubifex"),
    ("Mollusca", "Gastropoda", "Hygrophila", "Lymnaeidae", "Radix", "Radix balthica"),
]
INVALID = ["Baetis sp.", "Baetis sp. 1", "Chironomidae gen. sp.", "Gammarus cf. pulex",
           "Diptera sp. BIOUG1234-A01", "Simulium (Nevermannia) sp.", "None", ""]
PCTS = [100.0, 99.5, 98.0, 97.0, 96.99, 96.0, 95.0, 94.99, 92.0, 90.0, 89.9, 87.0,
        85.0, 84.99, 80.0, 75.0, 74.9, 60.0, 50.0, 49.0]
METHODS = [None, "Morphology", "BOLD ID Engine", "BIN Taxonomy Match", "Tree based", "DNA barcoding",
           "morphological", "BIN"]


def synthetic_tables(n, seed):
    rng = random.Random(seed)
    for s in range(n):
        sid = f"seq{s + 1}"
        if rng.random() < 0.04:
            rows = [dict(id=sid, **{l: "no-match" for l in LEVELS}, pct_identity=0.0, status="",
                         bin_uri="", process_id="")]
        else:
            k = rng.randint(1, 30)
            focus = rng.sample(TREE, rng.randint(1, 3))
            base = rng.choice(PCTS[: rng.randint(3, len(PCTS))])
            rows = []
            for _ in range(k):
                taxon = list(rng.choice(focus))
                # missing or invalid names at random levels (deeper levels more often)
                for j in range(6):
                    u = rng.random()
                    if u < 0.03 * (j + 1):
                        taxon[j] = None
                    elif u < 0.05 * (j + 1) and j >= 4:
                        taxon[j] = rng.choice(INVALID)
                pct = rng.choice([base, base, rng.choice(PCTS)])
                rows.append(dict(id=sid, **dict(zip(LEVELS, taxon)), pct_identity=pct,
                                 status=rng.choice(["public", "public", "private"]),
                                 bin_uri=rng.choice([None, "BOLD:AAA0001", "BOLD:AAB1234", "BOLD:ACX9999"]),
                                 process_id=f"P{rng.randint(1, 50)}",
                                 identification_method=rng.choice(METHODS)))
            # DuckDB orders by pct_identity DESC; ties come in some order: shuffle, then a stable sort
            rng.shuffle(rows)
            rows.sort(key=lambda r: -r["pct_identity"])
        df = pd.DataFrame(rows)
        for c in METADATA:
            if c not in df.columns:
                df[c] = None
        df["processid"] = df.get("process_id")
        df["coord"] = None
        df["fasta_order"] = s
        yield clean_dataframe(df)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "real":
        run(list(real_tables(sys.argv[2], sys.argv[3])), Path(sys.argv[4]))
    elif mode == "metadata":
        table = pd.read_csv(sys.argv[3], sep="\t", dtype=str, keep_default_na=False)
        ids = sorted(set(table["process_id"]) - {""})
        with duckdb.connect(sys.argv[2], read_only=True) as con:
            meta = con.execute("SELECT processid, identification_method FROM bold_public "
                               "WHERE processid IN (SELECT UNNEST(?))", [ids]).df()
        write(meta.to_dict("records"), ["processid", "identification_method"], Path(sys.argv[4]))
        print(f"{len(ids)} process IDs, {len(meta)} metadata rows -> {sys.argv[4]}")
    elif mode == "synthetic":
        n = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
        seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
        run(list(synthetic_tables(n, seed)), Path(sys.argv[2]))
    else:
        sys.exit(__doc__)
