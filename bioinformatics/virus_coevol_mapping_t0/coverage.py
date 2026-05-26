import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 15,
    "axes.titlesize": 15,
    "axes.labelsize": 15,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "axes.linewidth": 1.2,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})

# ─── INPUT ─────────────────────────────────────────────────────────────

coverage_file = "bioinformatics/virus_coevol_mapping_t0/virus_coverage_mean_1000.tsv"

# ─── LOAD ──────────────────────────────────────────────────────────────

df = pd.read_csv(
    coverage_file,
    sep="\t",
    header=None
)

df.columns = ["chrom", "start", "end", "coverage"]

# Replace missing values if any
df["coverage"] = pd.to_numeric(df["coverage"], errors="coerce")
df["coverage"] = df["coverage"].fillna(0.1)

# Avoid log(0)
df.loc[df["coverage"] <= 0, "coverage"] = 0.1

# ─── POSITION IN KBP ──────────────────────────────────────────────────

df["mid"] = (df["start"] + df["end"]) / 2
df["mid_kbp"] = df["mid"] / 1000

genome_size_kbp = df["end"].max() / 1000

# ─── FIGURE ───────────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(15, 4.5))

# Scatter plot
ax.scatter(
    df["mid_kbp"],
    df["coverage"],
    s=12,
    linewidths=0,
    rasterized=True,
)

# ─── X AXIS ────────────────────────────────────────────────────────────

# Major ticks every 20 kb
major_xticks = np.arange(
    0,
    genome_size_kbp + 20,
    20
)

ax.set_xticks(major_xticks)

ax.set_xlim(0, genome_size_kbp)

ax.tick_params(axis="x", which="major", length=6, width=1.0)
ax.tick_params(axis="x", which="minor", length=3, width=0.8)

# ─── LABELS ────────────────────────────────────────────────────────────

ax.set_xlabel("Genome position (kbp)")
ax.set_ylabel("Mean coverage per 1000 bp")

# ─── STYLE ─────────────────────────────────────────────────────────────

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.grid(
    axis="y",
    which="major",
    linestyle="--",
    linewidth=0.5,
    alpha=0.6
)

plt.tight_layout()

# ─── SAVE ──────────────────────────────────────────────────────────────

plt.savefig(
    "bioinformatics/virus_coevol_mapping_t0/virus_coverage_1000_log.png",
    dpi=300
)

plt.savefig(
    "bioinformatics/virus_coevol_mapping_t0/virus_coverage_1000_log.pdf"
)

plt.show()

print("Plots saved")