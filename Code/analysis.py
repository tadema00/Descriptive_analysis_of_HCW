"""Lassa fever surveillance analysis, Nigeria 2018 to epidemiological week 40 of 2026.

Reproduces the numbers in the manuscript:
  * Table 1/2/6 ratios and pooled 2020-2025 values (computed from combined counts)
  * Table 7 matched weeks 1-40 comparison (Fisher exact test, risk ratio with
    log-method 95% CI, risk difference with Wald 95% CI)
  * Supplementary Table S2 Cochran-Armitage trend tests over complete years
Run:  python analysis.py   (writes CSVs to outputs/)
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).parent
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

df = pd.read_csv(ROOT / "data" / "annual_indicators.csv").set_index("year")
matched = pd.read_csv(ROOT / "data" / "matched_weeks_1_40.csv").set_index("year")

# Bonferroni guide: nine reported tests (2 matched-period + 7 trend) as one family
N_TESTS = 9
ALPHA = 0.05 / N_TESTS  # 0.0056

# ---- indicator definitions: (numerator, denominator) as functions of a row ------
INDICATORS = {
    "CFR (deaths / confirmed)": lambda r: (r.deaths, r.confirmed),
    "Confirmed / suspected": lambda r: (r.confirmed, r.suspected),
    "HCW per confirmed (descriptive only)": lambda r: (r.hcw, r.confirmed),
    "Reported management ratio, managed / (confirmed + probable)": lambda r: (r.managed, r.confirmed + r.probable),
    "Follow-up completed / listed": lambda r: (r.completed, r.contacts_listed),
    "Lost to follow-up / listed": lambda r: (r.lost, r.contacts_listed),
    "Symptomatic / listed": lambda r: (r.symptomatic, r.contacts_listed),
    "Positive / listed": lambda r: (r.positive, r.contacts_listed),
}
COMPLETE = list(range(2018, 2026))          # complete years used for trend tests
POOLED = list(range(2020, 2026))            # pooled comparison years
TREND_YEARS = {"Reported management ratio, managed / (confirmed + probable)": list(range(2020, 2026))}


def ratio(a, n):
    return 100 * a / n


def pooled(f, years):
    num = den = 0
    for y in years:
        a, n = f(df.loc[y])
        num += a
        den += n
    return num, den  # combined counts, not a mean of annual percentages


# ---- 1. descriptive table (Table 6) ----------------------------------------------
rows = []
for name, f in INDICATORS.items():
    a26, n26 = f(df.loc[2026])
    a25, n25 = f(df.loc[2025])
    ap, npool = pooled(f, POOLED)
    rows.append({
        "indicator": name,
        "2026 (%)": ratio(a26, n26), "2025 (%)": ratio(a25, n25), "pooled 2020-25 (%)": ratio(ap, npool),
        "diff vs 2025 (points)": ratio(a26, n26) - ratio(a25, n25),
        "diff vs pooled (points)": ratio(a26, n26) - ratio(ap, npool),
    })
desc = pd.DataFrame(rows).round(2)
desc.to_csv(OUT / "table6_descriptive.csv", index=False)

# annual ratios for every year (Tables 1, 2; Figures 2-4)
annual = pd.DataFrame({name: [ratio(*f(df.loc[y])) if not pd.isna(f(df.loc[y])[0]) else np.nan
                              for y in df.index] for name, f in INDICATORS.items()}, index=df.index).round(2)
annual["contacts per confirmed case"] = (df.contacts_listed / df.confirmed).round(1)
annual.to_csv(OUT / "annual_ratios.csv")


# ---- 2. matched-period comparison (Table 7) --------------------------------------
def compare(a1, n1, a0, n0):
    """2026 (a1/n1) versus 2025 (a0/n0): Fisher p, RR (log CI), risk difference (Wald CI)."""
    p1, p0 = a1 / n1, a0 / n0
    rr = p1 / p0
    se_log = np.sqrt(1 / a1 - 1 / n1 + 1 / a0 - 1 / n0)
    rd = p1 - p0
    se_rd = np.sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
    _, p = stats.fisher_exact([[a1, n1 - a1], [a0, n0 - a0]])  # two-sided
    return {"2026 (%)": 100 * p1, "2025 (%)": 100 * p0,
            "RR": rr, "RR 95% CI low": rr * np.exp(-1.96 * se_log), "RR 95% CI high": rr * np.exp(1.96 * se_log),
            "RD (points)": 100 * rd, "RD 95% CI low": 100 * (rd - 1.96 * se_rd), "RD 95% CI high": 100 * (rd + 1.96 * se_rd),
            "Fisher p": p, "passes Bonferroni guide": p < ALPHA}


m26, m25 = matched.loc[2026], matched.loc[2025]
t7 = pd.DataFrame({
    "CFR (deaths / confirmed)": compare(m26.deaths, m26.confirmed, m25.deaths, m25.confirmed),
    "Confirmed / suspected": compare(m26.confirmed, m26.suspected, m25.confirmed, m25.suspected),
}).T
t7.round(4).to_csv(OUT / "table7_matched_weeks_1_40.csv")
counts = pd.DataFrame({k: {"2026": m26[k], "2025": m25[k], "difference": m26[k] - m25[k],
                           "relative change (%)": round(100 * (m26[k] / m25[k] - 1), 1)}
                       for k in ["suspected", "confirmed", "deaths", "hcw_weekly_notifications"]}).T
counts.to_csv(OUT / "table7_counts.csv")


# ---- 3. Cochran-Armitage trend tests (Supplementary Table S2) ----------------------
def cochran_armitage(events, totals, scores):
    """Two-sided Cochran-Armitage test for linear trend in proportions (score = year)."""
    r, n, x = (np.asarray(v, float) for v in (events, totals, scores))
    pbar = r.sum() / n.sum()
    num = (x * (r - n * pbar)).sum()
    den = np.sqrt(pbar * (1 - pbar) * ((n * x**2).sum() - (n * x).sum() ** 2 / n.sum()))
    z = num / den
    return z, 2 * stats.norm.sf(abs(z))


trend_rows = []
for name, f in INDICATORS.items():
    if name.startswith("HCW"):      # not tested: HCW infections are not necessarily a subset of confirmed cases
        continue
    yrs = TREND_YEARS.get(name, COMPLETE)
    pairs = [f(df.loc[y]) for y in yrs]
    z, p = cochran_armitage([a for a, _ in pairs], [n for _, n in pairs], yrs)
    trend_rows.append({"indicator": name, "years": f"{yrs[0]}-{yrs[-1]}", "z": round(z, 2), "p": p,
                       "direction": "falling" if z < -1.96 else "rising" if z > 1.96 else "no trend",
                       "passes Bonferroni guide": p < ALPHA})
pd.DataFrame(trend_rows).to_csv(OUT / "tableS2_trend_tests.csv", index=False)

# ---- 4. state concentration (Table 3) --------------------------------------------
states = pd.DataFrame({"state": ["Ondo", "Bauchi", "Taraba", "Edo", "Benue"],
                       "confirmed": [406, 265, 136, 122, 62], "deaths": [95, 36, 43, 27, 16],
                       "suspected": [2500, 1206, 391, 1941, 454]})
total_conf = df.loc[2026, "confirmed"]
states["share of national confirmed (%)"] = (100 * states.confirmed / total_conf).round(1)
states["CFR (%)"] = (100 * states.deaths / states.confirmed).round(1)
states["confirmed / suspected (%)"] = (100 * states.confirmed / states.suspected).round(1)
states.to_csv(OUT / "table3_leading_states.csv", index=False)
print(f"Five leading states: {states.confirmed.sum()} of {total_conf} = {100*states.confirmed.sum()/total_conf:.1f}%")

# ---- console summary --------------------------------------------------------------
pd.set_option("display.width", 200, "display.max_columns", 20)
print(f"Bonferroni guide: 0.05/{N_TESTS} = {ALPHA:.4f}\n")
print("Table 6 (descriptive):\n", desc.to_string(index=False), "\n")
print("Table 7 (matched weeks 1-40):\n", t7.round(4).to_string(), "\n")
print("Table 7 counts:\n", counts.to_string(), "\n")
print("Supplementary S2 (trend tests):\n", pd.DataFrame(trend_rows).to_string(index=False))
print("Totals:", df.suspected.sum(), df.confirmed.sum(), df.probable.sum(), df.deaths.sum())
