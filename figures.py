"""Figures 1-4 of the manuscript. Run after analysis.py:  python figures.py"""
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).parent
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
df = pd.read_csv(ROOT / "data" / "annual_indicators.csv").set_index("year")
GREY, BLUE, ORANGE, NAVY = "#c4c4b8", "#2878d4", "#d95f02", "#000080"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
labels = lambda ys: [f"{y}\n(wk 40)" if y == 2026 else str(y) for y in ys]


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=200)
    plt.close(fig)


# Figure 1: suspected and confirmed cases (2026 is a partial year: hollow marker, dashed segment)
fig, ax = plt.subplots(figsize=(8, 4.5))
for col, c, lab in [("suspected", NAVY, "Suspected"), ("confirmed", ORANGE, "Confirmed")]:
    ax.plot(df.index[:-1], df[col][:-1], "-o", color=c, label=lab)
    ax.plot(df.index[-2:], df[col][-2:], "--", color=c)
    ax.plot(df.index[-1], df[col].iloc[-1], "o", mfc="white", color=c)
    for y, v in df[col].items():
        ax.annotate(f"{v:,}", (y, v), textcoords="offset points", xytext=(0, 7), ha="center", color=c, fontsize=9)
ax.set_xticks(df.index); ax.set_xticklabels(labels(df.index))
ax.set_ylabel("Cases"); ax.legend(frameon=False)
ax.set_title("2026 is a partial year (to epidemiological week 40)", fontsize=10, loc="left")
save(fig, "figure1_cases.png")

# Figure 2: reported treatment-centre management ratio (managed / (confirmed + probable))
d = df.loc[2020:]
ratio = 100 * d.managed / (d.confirmed + d.probable)
fig, ax = plt.subplots(figsize=(8, 4.5))
bars = ax.bar(labels(d.index), ratio, color=[GREY] * (len(d) - 1) + [BLUE], width=0.65)
for b, v, m in zip(bars, ratio, d.managed):
    ax.text(b.get_x() + b.get_width() / 2, v + 1, f"{v:.1f}%\n({m:,} managed)", ha="center", fontsize=9)
ax.set_ylim(0, 115); ax.set_ylabel("Managed / (confirmed + probable), %")
ax.set_title("Reported treatment-centre management ratio; 2026 cumulative to week 40", fontsize=10, loc="left")
save(fig, "figure2_management_ratio.png")

# Figure 3: HCW infections and their descriptive ratio per confirmed case
d = df.loc[2020:]
share = 100 * d.hcw / d.confirmed
fig, ax = plt.subplots(figsize=(8, 4.5))
bars = ax.bar(range(len(d)), d.hcw, color=[GREY] * (len(d) - 1) + [BLUE], width=0.65)
for i, (b, v, s) in enumerate(zip(bars, d.hcw, share)):
    ax.text(i, v + 1, str(v), ha="center")
ax.set_xticks(range(len(d)))
ax.set_xticklabels([f"{l}\n{s:.1f}%" for l, s in zip(labels(d.index), share)])
ax.set_ylabel("HCW infections")
ax.set_title("HCW infections; % = descriptive ratio per 100 confirmed cases (not a risk). 2025 full year; 2026 to week 40",
             fontsize=9, loc="left")
save(fig, "figure3_hcw.png")

# Figure 4: symptomatic and positive contacts as a share of ALL contacts listed (not contacts tested)
d = df.loc[2020:]
fig, ax = plt.subplots(figsize=(8, 4.5))
w = 0.38; x = range(len(d))
s = 100 * d.symptomatic / d.contacts_listed; p = 100 * d.positive / d.contacts_listed
ax.bar([i - w / 2 for i in x], s, w, color=BLUE, label="Symptomatic contacts")
ax.bar([i + w / 2 for i in x], p, w, color=GREY, label="Positive contacts")
for i, (a, b) in enumerate(zip(s, p)):
    ax.text(i - w / 2, a + 0.1, f"{a:.1f}", ha="center", fontsize=9); ax.text(i + w / 2, b + 0.1, f"{b:.1f}", ha="center", fontsize=9)
ax.set_xticks(list(x)); ax.set_xticklabels([f"{l}\n{n:,} listed" for l, n in zip(labels(d.index), d.contacts_listed)])
ax.set_ylabel("% of contacts listed"); ax.legend(frameon=False, loc="upper left")
ax.set_title("Denominator is all contacts listed, not contacts tested; 2026 to week 40", fontsize=10, loc="center")
save(fig, "figure4_contacts.png")
print("Figures written to", OUT)
