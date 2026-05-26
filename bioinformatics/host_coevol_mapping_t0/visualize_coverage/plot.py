import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import numpy as np

matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 15,
    "axes.titlesize": 15,
    "axes.labelsize": 15,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "axes.linewidth": 1.2,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})

# ─── Chromosome ordering & naming ───────────────────────────────────────────
chrom_order = [
    "CP001574.1",  # Chr 1
    "CP001323.1",  # Chr 2
    "CP001324.1",  # Chr 3
    "CP001325.1",  # Chr 4
    "CP001326.1",  # Chr 5
    "CP001327.1",  # Chr 6
    "CP001328.1",  # Chr 7
    "CP001575.1",  # Chr 8
    "CP001329.1",  # Chr 9
    "CP001576.1",  # Chr 10
    "CP001330.1",  # Chr 11
    "CP001577.1",  # Chr 12
    "CP001331.1",  # Chr 13
    "CP001332.1",  # Chr 14
    "CP001333.1",  # Chr 15
    "CP001334.1",  # Chr 16
    "CP001335.1",  # Chr 17
]

chrom_name_map = {
    "CP001574.1": "Chr 1",
    "CP001323.1": "Chr 2",
    "CP001324.1": "Chr 3",
    "CP001325.1": "Chr 4",
    "CP001326.1": "Chr 5",
    "CP001327.1": "Chr 6",
    "CP001328.1": "Chr 7",
    "CP001575.1": "Chr 8",
    "CP001329.1": "Chr 9",
    "CP001576.1": "Chr 10",
    "CP001330.1": "Chr 11",
    "CP001577.1": "Chr 12",
    "CP001331.1": "Chr 13",
    "CP001332.1": "Chr 14",
    "CP001333.1": "Chr 15",
    "CP001334.1": "Chr 16",
    "CP001335.1": "Chr 17",
}

# Alternating two-tone palette (only nuclear chromosomes now)
palette = [
    "#2166ac", "#4393c3",
    "#d6604d", "#f4a582",
    "#1a9641", "#74c476",
    "#762a83", "#af8dc3",
    "#e08214", "#fec44f",
    "#35978f", "#99d8c9",
    "#543005", "#a6611a",
    "#de77ae", "#f1b6da",
]

chrom_color_map = {chrom: palette[i % len(palette)] for i, chrom in enumerate(chrom_order)}

# ─── Paths ───────────────────────────────────────────────────────────────────
path_input  = "bioinformatics/host_coevol_mapping_t0/visualize_coverage/genotoul_output"
path_output = "bioinformatics/host_coevol_mapping_t0/visualize_coverage/python_output"

# ─── Load data ───────────────────────────────────────────────────────────────
df = pd.read_csv(
    f"{path_input}/Mc_cat_mem2_sorted_coverage1000_sorted",
    sep="\t",
    header=None
)

df.columns = ["chrom", "start", "end", "value", "col4", "length", "ratio"]

# ─── REMOVE organelles (mitochondrie + chloroplaste) ────────────────────────
excluded = {"FJ859351.1", "FJ858267.1"}
df = df[~df["chrom"].isin(excluded)].copy()

# ─── Sort chromosomes ────────────────────────────────────────────────────────
df["chrom"] = pd.Categorical(df["chrom"], categories=chrom_order, ordered=True)
df = df.sort_values(["chrom", "start"]).reset_index(drop=True)

# ─── Compute genome-wide x positions ─────────────────────────────────────────
df["mid"] = (df["start"] + df["end"]) / 2

x_positions   = []
chrom_offsets = {}
offset        = 0

for chrom, sub in df.groupby("chrom", sort=False):
    chrom_offsets[chrom] = offset
    x_positions.extend(sub["mid"] + offset)
    offset += sub["end"].max()

df["x"] = x_positions

# ─── Chromosome boundary lines ───────────────────────────────────────────────
boundaries = []
for chrom, sub in df.groupby("chrom", sort=False):
    boundaries.append(chrom_offsets[chrom] + sub["end"].max())

# ─── Figure ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(16, 4.5))

for chrom, sub in df.groupby("chrom", sort=False):
    ax.scatter(
        sub["x"],
        sub["value"],
        s=3,
        color=chrom_color_map[chrom],
        linewidths=0,
        rasterized=True,
        label=chrom_name_map.get(chrom, chrom),
        zorder=2,
    )

# Vertical separators
for b in boundaries[:-1]:
    ax.axvline(b, color="#cccccc", linewidth=0.6, zorder=1)

# ─── X axis labels ───────────────────────────────────────────────────────────
xticks = []
xticklabels = []

for chrom, sub in df.groupby("chrom", sort=False):
    center = chrom_offsets[chrom] + sub["end"].max() / 2
    xticks.append(center)
    xticklabels.append(chrom_name_map.get(chrom, chrom))

ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels, rotation=90, ha="center", va="top")

ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True)

# ─── Y axis ──────────────────────────────────────────────────────────────────
ax.yaxis.set_major_locator(ticker.AutoLocator())
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
ax.tick_params(axis="y", which="major", length=5, width=1.0)
ax.tick_params(axis="y", which="minor", length=2.5, width=0.7)

# ─── Labels & formatting ─────────────────────────────────────────────────────
ax.set_ylabel("Coverage (reads per 1000bp)", labelpad=8)
#ax.set_title("Genome-wide read coverage — window size 1,000 bp", pad=10, fontweight="bold")

padding = df.groupby("chrom")["end"].max().iloc[-1] * 0.6
ax.set_xlim(0, df["x"].max() + padding)

ax.set_yscale("log")
ax.spines[["top", "right"]].set_visible(False)
ax.spines["bottom"].set_linewidth(1.2)
ax.spines["left"].set_linewidth(1.2)

ax.grid(axis="y", which="major", linestyle="--", linewidth=0.5, alpha=0.6, zorder=0)

fig.subplots_adjust(bottom=0.18)

# ─── Save ─────────────────────────────────────────────────────────────────────
out_png = f"{path_output}/coverage_1000.png"
out_pdf = f"{path_output}/coverage_1000.pdf"

plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.show()

print("Saved:", out_png, "&", out_pdf)