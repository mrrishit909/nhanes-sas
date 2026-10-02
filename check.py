"""The project's check: recompute the SAS results in Python, and compare both with CDC's published figures.

    ./venv/bin/python check.py

1. Every survey-weighted proportion the SAS program reports (results/estimates.csv, parsed from the SAS log) is
   recomputed here from the same CDC files: the estimate, and its standard error by Taylor linearisation with
   strata and PSUs (the method PROC SURVEYMEANS uses).
2. The headline figures round to the ones in NCHS Data Brief No. 511.
3. The odds ratios from PROC SURVEYLOGISTIC match a weighted logistic regression fitted here (point estimates
   do not depend on the design; only their standard errors do).
"""
import ssl
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

HERE = Path(__file__).parent
RAW = HERE / "data" / "raw"
R = HERE / "results"
FILES = ["DEMO_L", "BPXO_L", "BPQ_L", "BMX_L", "HIQ_L"]
URL = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/{}.xpt"
# NCHS Data Brief No. 511 (October 2024), August 2021 - August 2023, percent
PUBLISHED = {("htn", "All"): 47.7, ("htn", "Men"): 50.8, ("htn", "Women"): 44.6, ("htn", "18-39"): 23.4, ("htn", "40-59"): 52.5,
             ("htn", "60+"): 71.6, ("aware", "All"): 59.2, ("aware", "Men"): 55.2, ("aware", "Women"): 63.6,
             ("aware", "18-39"): 27.2, ("aware", "40-59"): 56.7, ("aware", "60+"): 73.7, ("meds", "All"): 51.2,
             ("meds", "Men"): 46.7, ("meds", "Women"): 56.1, ("meds", "18-39"): 13.9, ("meds", "40-59"): 47.1,
             ("meds", "60+"): 69.1, ("controlled", "All"): 20.7, ("controlled", "Men"): 18.9, ("controlled", "Women"): 22.8,
             ("controlled", "18-39"): 4.5, ("controlled", "40-59"): 18.1, ("controlled", "60+"): 29.2}
PUBLISHED_AGEADJ = {"All": 44.5, "Men": 48.8, "Women": 40.1}


def load():
    RAW.mkdir(parents=True, exist_ok=True)
    ctx = ssl.create_default_context(cafile="/etc/ssl/cert.pem")
    for f in FILES:
        if not (RAW / f"{f}.xpt").exists():
            req = urllib.request.Request(URL.format(f), headers={"User-Agent": "Mozilla/5.0"})
            (RAW / f"{f}.xpt").write_bytes(urllib.request.urlopen(req, context=ctx, timeout=120).read())
    d = pd.read_sas(RAW / "DEMO_L.xpt", format="xport")
    for f in FILES[1:]:
        d = d.merge(pd.read_sas(RAW / f"{f}.xpt", format="xport"), on="SEQN", how="left")
    d["sbp"] = d[["BPXOSY1", "BPXOSY2", "BPXOSY3"]].mean(axis=1)
    d["dbp"] = d[["BPXODI1", "BPXODI2", "BPXODI3"]].mean(axis=1)
    d["adult"] = (d.RIDAGEYR >= 18) & d.sbp.notna() & d.dbp.notna() & (d.WTMEC2YR > 0) & d.RIDEXPRG.ne(1)
    d["meds"] = d.BPQ150.eq(1).astype(float)
    d["htn"] = ((d.sbp >= 130) | (d.dbp >= 80) | d.meds.eq(1)).astype(float)
    d["aware"] = np.where(d.BPQ020.isin([1, 2]), d.BPQ020.eq(1), np.nan)
    d["controlled"] = ((d.sbp < 130) & (d.dbp < 80)).astype(float)
    d["hyp"] = d.adult & d.htn.eq(1)
    d["treated"] = d.hyp & d.meds.eq(1)
    d["unaware"] = np.where(d.hyp & d.aware.notna(), 1 - d.aware, np.nan)
    d["uncontrolled"] = np.where(d.treated, 1 - d.controlled, np.nan)
    d["agegrp"] = np.select([d.RIDAGEYR < 40, d.RIDAGEYR < 60], ["18-39", "40-59"], "60+")
    d["sex"] = np.where(d.RIAGENDR == 1, "Men", "Women")
    d["race"] = d.RIDRETH3.map({1: "Hispanic", 2: "Hispanic", 3: "White", 4: "Black", 6: "Asian"}).fillna("Other")
    p = d.INDFMPIR
    d["income"] = np.where(p.isna(), None, np.select([p < 2, p < 4], ["<2x", "2-4x"], "4x+"))      # unknowns left out, as in SAS
    d["insured"] = d.HIQ011.map({1: "yes", 2: "no"})
    b = d.BMXBMI
    d["bmi"] = np.where(b.isna(), None, np.select([b < 25, b < 30], ["<25", "25-30"], "30+"))
    return d


def domain_mean(d, y, in_domain):
    """Ratio estimate of a domain mean and its Taylor-linearised SE (with-replacement PSUs within strata).
    Like SURVEYMEANS with NOMCAR, cases in the domain with a missing y are treated as a further domain."""
    dom = in_domain & d[y].notna()
    w = d.WTMEC2YR.where(dom, 0.0)
    yv = d[y].fillna(0.0)
    est = (w * yv).sum() / w.sum()
    z = w * (yv - est) / w.sum()
    t = z.groupby([d.SDMVSTRA, d.SDMVPSU]).sum()
    var = 0.0
    for _, th in t.groupby(level=0):
        n = len(th)
        var += n / (n - 1) * ((th - th.mean()) ** 2).sum()
    return est, np.sqrt(var), int(dom.sum())


def main():
    d = load()
    sas = pd.read_csv(R / "estimates.csv")
    domains = {"htn": "adult", "aware": "hyp", "meds": "hyp", "controlled": "hyp"}
    n_checked = 0
    for _, row in sas.iterrows():
        dom = d[domains[row["measure"]]]
        if row["level"] != "All":
            for part in row["level"].split(" "):
                dom = dom & (d.sex.eq(part) | d.agegrp.eq(part))
        est, se, n = domain_mean(d, row["measure"], dom)
        assert n == row["n"], (row.to_dict(), n)
        assert abs(est - row["mean"]) < 1e-6 and abs(se - row["stderr"]) < 1e-6, (row.to_dict(), est, se)
        if (row["measure"], row["level"]) in PUBLISHED:
            assert round(100 * row["mean"], 1) == PUBLISHED[(row["measure"], row["level"])], row.to_dict()
        n_checked += 1
    adj = pd.read_csv(R / "ageadj.csv").set_index("group")["adj_mean"]
    for g, v in PUBLISHED_AGEADJ.items():
        assert round(100 * adj[g], 1) == v, (g, adj[g])

    odds = pd.read_csv(R / "odds.csv")
    ref = {"agegrp": "18-39", "sex": "Women", "race": "White", "income": "4x+", "insured": "yes", "bmi": "<25"}   # 12 parameters
    for y, dom in [("unaware", "hyp"), ("uncontrolled", "treated")]:
        m = d[d[dom] & d[y].notna() & d[list(ref)].notna().all(axis=1)]
        X = pd.concat([pd.get_dummies(m[v], prefix=v).drop(columns=f"{v}_{r}") for v, r in ref.items()], axis=1).astype(float)
        fit = sm.GLM(m[y], sm.add_constant(X), family=sm.families.Binomial(), var_weights=m.WTMEC2YR).fit()
        for _, row in odds[odds["model"] == y].iterrows():
            var, level = row["effect"].split(" ", 1)
            level = level.split(" vs ")[0]
            assert abs(np.exp(fit.params[f"{var}_{level}"]) - row["oddsratioest"]) < 1e-3, (y, row.to_dict())
            n_checked += 1
    print(f"OK: {n_checked} SAS results recomputed in Python (estimates and design-based SEs to 1e-6, odds ratios to 1e-3); "
          f"{len(PUBLISHED)} figures and 3 age-adjusted rates match NCHS Data Brief 511")


if __name__ == "__main__":
    main()
