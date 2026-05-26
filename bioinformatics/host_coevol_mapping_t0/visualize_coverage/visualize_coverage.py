import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

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
    "FJ859351.1",  # Mitochondrie
    "FJ858267.1"   # Chloroplaste
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
    "FJ859351.1": "Mitochondrie",
    "FJ858267.1": "Chloroplaste"
}

colors = cm.tab20.colors
chrom_color_map = {chrom: colors[i % len(colors)] for i, chrom in enumerate(chrom_order)}

path_input = "bioinformatics/host_coevol_mapping_t0/visualize_coverage/genotoul_output"
path_output = "bioinformatics/host_coevol_mapping_t0/visualize_coverage/python_output"

# =========================
# Lecture des données
# =========================
df1000 = pd.read_csv(f"{path_input}/Sputnik1_coverage1000", sep="\t", header=None)
df1000.columns = ["chrom", "start", "end", "value", "col4", "length", "ratio"]

df10000 = pd.read_csv(f"{path_input}/Sputnik1_coverage1000", sep="\t", header=None)
df10000.columns = ["chrom", "start", "end", "value", "col4", "length", "ratio"]

# =========================
# Tri des chromosomes
# =========================
df1000["chrom"] = pd.Categorical(df1000["chrom"], categories=chrom_order, ordered=True)
df1000 = df1000.sort_values(by=["chrom", "start"]).reset_index(drop=True)

df10000["chrom"] = pd.Categorical(df10000["chrom"], categories=chrom_order, ordered=True)
df10000 = df10000.sort_values(by=["chrom", "start"]).reset_index(drop=True)

# =========================
# Position des bins
# =========================
df1000["mid"] = (df1000["start"] + df1000["end"]) / 2

x_positions = []
offset = 0
chrom_offsets = {}

for chrom, subdf1000 in df1000.groupby("chrom", sort=False):
    chrom_offsets[chrom] = offset
    x_positions.extend(subdf1000["mid"] + offset)
    offset += subdf1000["end"].max()

df1000["x"] = x_positions


df10000["mid"] = (df10000["start"] + df10000["end"]) / 2

x_positions = []
offset = 0
chrom_offsets = {}

for chrom, subdf10000 in df10000.groupby("chrom", sort=False):
    chrom_offsets[chrom] = offset
    x_positions.extend(subdf10000["mid"] + offset)
    offset += subdf10000["end"].max()

df10000["x"] = x_positions

# =========================
# 7. Plot
# =========================

plt.subplot(2,1,1)

for chrom, subdf1000 in df1000.groupby("chrom", sort=False):
    plt.plot(
        subdf1000["x"],
        subdf1000["value"],
        marker="o",
        linestyle="",
        markersize=0.5,
        color=chrom_color_map[chrom],
        label=chrom_name_map.get(chrom, chrom)
    )

# Axe X : centres des chromosomes
xticks = []
xticklabels = []

for chrom, subdf1000 in df1000.groupby("chrom", sort=False):
    center = chrom_offsets[chrom] + subdf1000["end"].max() / 2
    xticks.append(center)
    xticklabels.append(chrom_name_map.get(chrom, chrom))

plt.xticks(xticks, xticklabels, rotation=45)
plt.ylabel("Coverage")
plt.title("Window size : 1000")



plt.subplot(2,1,2)

for chrom, subdf10000 in df10000.groupby("chrom", sort=False):
    plt.plot(
        subdf10000["x"],
        subdf10000["value"],
        marker="o",
        linestyle="",
        markersize=0.5,
        color=chrom_color_map[chrom],
        label=chrom_name_map.get(chrom, chrom)
    )

# Axe X : centres des chromosomes
xticks = []
xticklabels = []

for chrom, subdf10000 in df10000.groupby("chrom", sort=False):
    center = chrom_offsets[chrom] + subdf10000["end"].max() / 2
    xticks.append(center)
    xticklabels.append(chrom_name_map.get(chrom, chrom))

plt.xticks(xticks, xticklabels, rotation=45)
plt.ylabel("Coverage")
plt.title("Window size : 10000")


plt.tight_layout()
plt.savefig(f"{path_output}/coverage_1000_10000.png", dpi=300, bbox_inches="tight")
plt.show()
