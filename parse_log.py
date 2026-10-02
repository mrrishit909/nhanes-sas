"""Turn the CSV| lines of a saved SAS log into results files.

    ./venv/bin/python parse_log.py sas/01_hypertension.log   -> results/estimates.csv, results/ageadj.csv, results/odds.csv
"""
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
COLS = {"estimate": ["measure", "level", "n", "mean", "stderr", "lowerclmean", "upperclmean"],
        "ageadj": ["group", "adj_mean", "adj_se"], "odds": ["model", "effect", "oddsratioest", "lowercl", "uppercl"]}
OUT = {"estimate": "estimates.csv", "ageadj": "ageadj.csv", "odds": "odds.csv"}


def main(log):
    rows = {k: [] for k in COLS}
    for line in Path(log).read_text().splitlines():
        if line.startswith("CSV|"):
            parts = line.strip().split("|")
            rows[parts[1]].append(parts[2:])
    (HERE / "results").mkdir(exist_ok=True)
    for k, r in rows.items():
        assert r, f"no {k} lines in the log"
        pd.DataFrame(r, columns=COLS[k]).to_csv(HERE / "results" / OUT[k], index=False)
        print(OUT[k], len(r), "rows")


if __name__ == "__main__":
    main(sys.argv[1])
