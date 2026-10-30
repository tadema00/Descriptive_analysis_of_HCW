# Lassa fever in Nigeria: surveillance indicators, 2018 to epidemiological week 40 of 2026

Code and data for the manuscript *Epidemiological Trends and Changes in Lassa Fever Surveillance Indicators in
Nigeria, 2018–2026: A Retrospective Descriptive Study*.

## Contents
| File | Purpose |
|---|---|
| `data/annual_indicators.csv` | Annual national counts from NCDC situation reports (week 52; week 53 in 2020; week 40 in 2026) |
| `data/matched_weeks_1_40.csv` | Cumulative weeks 1–40 values for 2025 and 2026 (HCW column = weekly notifications) |
| `analysis.py` | Descriptive ratios, matched-period tests, Cochran–Armitage trend tests, state table |
| `figures.py` | Figures 1–4 |
| `outputs/` | Generated tables (CSV) and figures (PNG) |

## Run
```bash
pip install -r requirements.txt
python analysis.py
python figures.py
```
Tested with Python 3 and SciPy 1.18.

## Methods in brief
* Pooled 2020–2025 values use combined counts, not means of annual percentages.
* Matched weeks 1–40 (2026 vs 2025): two-sided Fisher exact test, risk ratio (log-method 95% CI), risk difference (Wald 95% CI).
* Trend tests: Cochran–Armitage with year as the score, over complete years (2018–2025; 2020–2025 for management).
* Bonferroni guide: nine tests, α = 0.05/9 = 0.0056. This is exploratory, not prespecified.
* HCW infections per confirmed case is a descriptive ratio and is not tested (HCW infections are not necessarily a subset of confirmed cases).
* "Management ratio" is managed ÷ (confirmed + probable); the NCDC operational definition is not published.

## Data notes
Where sources disagreed, the NCDC SitRep value was used (see manuscript Section 2.5). Follow-up categories do not
add up to contacts listed in 2018–2024. 2026 values are provisional.
Source: Nigeria Centre for Disease Control and Prevention, Lassa fever situation reports (https://ncdc.gov.ng/diseases/sitreps).
