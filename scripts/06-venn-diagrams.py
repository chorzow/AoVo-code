"""
venn diagrams
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml
from matplotlib_venn import venn2, venn2_circles


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent

with open(ROOT_DIR / "config.yml") as f:
    CONFIG = yaml.safe_load(f)["default"]

DATA_OUT = ROOT_DIR / CONFIG["output_folder"] / "data"
IMAGES_OUT = ROOT_DIR / CONFIG["images_output_folder"]


def families_present(df: pd.DataFrame, samples: list[str]) -> set[str]:
    """Families with a non-zero count in at least one of samples."""
    if not samples:
        return set()
    return set(df.index[(df[samples] > 0).any(axis=1)])


def jaccard_index(set1: set[str], set2: set[str]) -> float:
    return len(set1 & set2) / len(set1 | set2)


def draw_venn(ax, set1, set2, labels, colors, title):
    v = venn2((set1, set2), set_labels=labels, ax=ax)
    for patch_id, color in zip(("10", "01", "11"), (*colors, None)):
        patch = v.get_patch_by_id(patch_id)
        if patch is None:
            continue
        if color is not None:
            patch.set_color(color)
        patch.set_alpha(0.7)
    venn2_circles((set1, set2), linewidth=0.5, ax=ax)
    ax.set_title(title, fontweight="bold", loc="left")


data_family = pd.read_csv(DATA_OUT / "16S_FamilyAOVO.tsv", sep="\t", index_col=0)
data_family = data_family.drop(columns=["total"], errors="ignore")

groups = {
    "S1_bulk": [c for c in data_family.columns if "S1_bulk" in c],
    "S1_rAo": [c for c in data_family.columns if "S1_rAo" in c],
    "S2_bulk": [c for c in data_family.columns if "S2_bulk" in c],
    "S2_rVo": [c for c in data_family.columns if "S2_rVo" in c],
}
families = {group: families_present(data_family, samples) for group, samples in groups.items()}

jaccard_sites = jaccard_index(families["S1_bulk"], families["S2_bulk"])
jaccard_s1 = jaccard_index(families["S1_bulk"], families["S1_rAo"])
jaccard_s2 = jaccard_index(families["S2_bulk"], families["S2_rVo"])
print(f"Jaccard index (S1_bulk vs S2_bulk, sites):  {jaccard_sites:.3f}")
print(f"Jaccard index (S1_bulk vs S1_rAo):          {jaccard_s1:.3f}")
print(f"Jaccard index (S2_bulk vs S2_rVo):          {jaccard_s2:.3f}")

comparisons = [
    ("A", families["S1_bulk"], families["S2_bulk"], ("S1_bulk", "S2_bulk"), ("#FF9999", "#99FF99")),
    ("B", families["S1_bulk"], families["S1_rAo"], ("S1_bulk", "S1_rAo"), ("#E41A1C", "#377EB8")),
    ("C", families["S2_bulk"], families["S2_rVo"], ("S2_bulk", "S2_rVo"), ("#4DAF4A", "#984EA3")),
]

fig = plt.figure(figsize=(10, 10))
gs = fig.add_gridspec(2, 2)
axes = [fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]

for ax, (title, set1, set2, labels, colors) in zip(axes, comparisons):
    draw_venn(ax, set1, set2, labels, colors, title)

fig.tight_layout()

IMAGES_OUT.mkdir(parents=True, exist_ok=True)
fig.savefig(IMAGES_OUT / "Figure_Venn_diagrams_fixed.pdf")
