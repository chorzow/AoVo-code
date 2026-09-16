# AoVo

16S rRNA amplicon analysis of *Astragalus olchonensis* / *Vicia olchonensis*
rhizosphere and bulk-soil samples.

## Prerequisites

- **Python 3.11+**
- [uv](https://docs.astral.sh/uv/)
- **R**, either via:
  - [uvr](https://github.com/nbafrank/uvr) (recommended — pins the R
    version and package library via `uvr.toml`/`uvr.lock`), or
  - a plain R installation with `Rscript` on PATH and this project's R
    packages (see `uvr.toml` for the list) already installed.

With uvr, install its R version and packages once:

```bash
uvr r install 4.5.3
uvr sync
```

Install the Python environment:

```bash
uv sync
```

## Running the pipeline

```bash
uv run main.py
```

This runs `scripts/01-preload.R` through `scripts/07-indicator-species.R` in
order, then `06-venn-diagrams.py` via the active Python interpreter. R
scripts run through `uvr run` if `uvr.toml` is present and `uvr` is
installed, falling back to plain `Rscript` otherwise. Force one explicitly
with `R_RUNNER`:

```bash
R_RUNNER=uvr uv run main.py       # always use uvr
R_RUNNER=rscript uv run main.py   # always use plain Rscript
```

Each R script sources `01-preload.R` for shared setup, so individual steps
can also be run on their own, e.g.:

```bash
uvr run scripts/05-taxonomic-analysis.R
# or, without uvr:
Rscript scripts/05-taxonomic-analysis.R
```

Outputs are written to `output/images/` (PDFs) and `output/data/`
(TSV/CSV tables), as configured in `config.yml`.

## Pipeline stages

| Script | Produces |
|---|---|
| `01-preload.R` | Shared setup: libraries, config, OTU/taxonomy/metadata loading, helper functions — sourced by every other script |
| `02-rarecurve.R` | Rarefaction curve of observed OTUs per sample group |
| `03-alpha-diversity.R` | Chao1/Shannon/Simpson diversity, boxplots, pairwise Wilcoxon tests |
| `04-beta-diversity.R` | NMDS, PCoA, ANOSIM, PERMANOVA on Bray-Curtis distance |
| `05-taxonomic-analysis.R` | Dominant Phylum/Family composition, variability analysis, Kruskal-Wallis/Dunn test, and a microshades composition plot |
| `06-venn-diagrams.py` | Venn diagrams of shared Families between sample groups |
| `07-indicator-species.R` | IndVal indicator-species analysis (`multipatt`), boxplots and heatmaps of top indicators |

## Configuration

`config.yml` holds every path, sample-group definition, and analysis
parameter (rarefaction step, top-N taxa, p-value threshold, indicator
thresholds, etc.) used by the pipeline — edit it rather than the scripts to
change these.

## Data inputs

Expected under `data/` (paths configurable in `config.yml`):

- `otu_table.tsv` — OTU counts (OTUs x samples)
- `taxonomy.tsv` — taxonomic assignment per OTU
- `metadata.csv` — environmental/sample metadata; rows are matched against `otu_table.tsv`'s columns
  via its `Sample` column (e.g. `S1_bulk1`)
