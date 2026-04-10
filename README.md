# constant-ph-utils

Analysis and plotting toolkit for **GROMACS constant-pH molecular dynamics** simulations.  
Reads λ-coordinate XVG files, computes protonation statistics across replicas, and generates publication-ready figures.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
  - [Shell Completion](#shell-completion)
- [Input Data](#input-data)
- [Usage](#usage)
  - [Quick Start](#quick-start)
  - [Full CLI Reference](#full-cli-reference)
  - [Examples](#examples)
- [Output Files](#output-files)
- [Module Reference](#module-reference)
  - [analyses.py](#analysespy)
  - [plot.py](#plotpy)
- [Notes on Histidine Handling](#notes-on-histidine-handling)
- [Multi-Chain (Homomeric) Systems](#multi-chain-homomeric-systems)
- [Profiling](#profiling)

---

## Overview

Constant-pH MD simulations in GROMACS produce per-residue λ-coordinate time series (`cphmd-coord-*.xvg`).
This toolkit loads those files in parallel, classifies frames as protonated (λ < 0.2) or deprotonated (λ > 0.8), and provides:

- **Lambda histograms** — distribution of λ values per residue.
- **Protonation time series** — cumulative protonation fraction over time.
- **Convergence plots** — mean ± SE across replicas as a function of simulation time.
- **Protonation fraction bar charts** — publication-ready panels grouped by residue type, with error bars from replica statistics.
- **Single-residue convergence** — focused convergence view for residues of interest.

All plots support both a `Debug` mode (coordinate IDs in titles) and a `Publication` mode (residue names, clean styling).

---

## Repository Structure

```
constant-ph-utils/
├── analyses.py                  # Data loading, statistics, time series
├── plot.py                      # All plotting functions
├── run_analysis_template.py     # CLI entry point
└── README.md
```

---

## Requirements

- **Python ≥ 3.10** (uses `int | slice | tuple` union syntax)
- **NumPy**
- **pandas**
- **matplotlib ≥ 3.3** (uses `GridSpec.subgridspec` and `bar_label`)
- **argcomplete** (optional, for shell tab-completion of CLI flags)

All other imports (`dataclasses`, `concurrent.futures`, `argparse`, `cProfile`, etc.) are part of the Python standard library.

---

## Installation

No installation step is required. Clone or copy the three Python files into a directory and run directly:

```bash
git clone <repo-url> constant-ph-utils
cd constant-ph-utils
```

Install the Python dependencies if not already available:

```bash
pip install numpy pandas matplotlib argcomplete
```

### Shell Completion

To enable tab-completion of all `--flags` in your terminal, run the appropriate activation for your shell:

**Bash** (add to `~/.bashrc`):

```bash
eval "$(register-python-argcomplete run_analysis_template.py)"
```

**Zsh** (add to `~/.zshrc`):

```bash
autoload -U bashcompinit
bashcompinit
eval "$(register-python-argcomplete run_analysis_template.py)"
```

**Global activation** (enables completion for all argcomplete-enabled scripts):

```bash
activate-global-python-argcomplete
```

After sourcing your shell config (`source ~/.bashrc` or `source ~/.zshrc`), you can type:

```bash
python run_analysis_template.py --<TAB><TAB>
```

and get a list of all available flags.

---

## Input Data

The toolkit expects the following directory layout:

```
project/
├── lambdareference.dat                # Lambda reference table
├── MD1/
│   └── analysis/
│       ├── cphmd-coord-1.xvg
│       ├── cphmd-coord-2.xvg
│       └── ...
└── MD2/
    └── analysis/
        ├── cphmd-coord-1.xvg
        ├── cphmd-coord-2.xvg
        └── ...
```

**`lambdareference.dat`** — whitespace-separated table with at least these columns:

| Column           | Description                            |
|------------------|----------------------------------------|
| `resname`        | Residue name (e.g. `GLUT`, `HSPT`)    |
| `resid`          | Residue number                         |
| `chain`          | Chain identifier                       |
| `coordinateFile` | Corresponding `cphmd-coord-*.xvg` file |

**`cphmd-coord-*.xvg`** — GROMACS-style XVG files with two columns: time and λ value. Lines starting with `#` or `@` are skipped.

---

## Usage

### Quick Start

Minimal invocation with two replicas:

```bash
python run_analysis_template.py \
  --lambdaref-path project/ \
  --paths-md project/MD1/analysis project/MD2/analysis
```

This runs in `Publication` mode (the default), producing all plot types.

### Full CLI Reference

```
usage: run_analysis_template.py [-h]
    --lambdaref-path LAMBDAREF_PATH
    --paths-md PATHS_MD [PATHS_MD ...]
    [--run-type {Debug,Publication}]
    [--plot-type {Debug,Publication}]
    [--plot-rows PLOT_ROWS]
    [--plot-cols PLOT_COLS]
    [--no-single-letter]
    [--npz-output]
    [--dpi DPI]
    [--threads THREADS]
    [--xvg-rows XVG_ROWS]
    [--chains CHAINS [CHAINS ...]]
    [--res-ids RES_IDS [RES_IDS ...]]
    [--profile]
```

| Argument | Default | Description |
|---|---|---|
| `--lambdaref-path` | *required* | Directory containing `lambdareference.dat` |
| `--paths-md` | *required* | One or more paths to MD analysis directories |
| `--run-type` | `Publication` | `Debug` = single-replica plots only; `Publication` = adds convergence + fraction plots |
| `--plot-type` | `Debug` | `Debug` = coord IDs in titles; `Publication` = residue names |
| `--plot-rows` | `21` | Rows in the overview plot grids |
| `--plot-cols` | `5` | Columns in the overview plot grids |
| `--no-single-letter` | *(off)* | Use three-letter amino acid codes (default is one-letter) |
| `--npz-output` | `False` | Save protonation data as `.npz` files |
| `--dpi` | `300` | DPI for all saved figures |
| `--threads` | `8` | Parallel threads for XVG loading |
| `--xvg-rows` | `2000000` | Max rows to read per XVG file |
| `--chains` | `None` | Chain identifiers for homomeric systems (e.g. `A B`) |
| `--res-ids` | `None` | Residue IDs for single-residue convergence plots |
| `--profile` | `False` | Enable cProfile profiling |

### Examples

**Debug mode — single replica, quick check:**

```bash
python run_analysis_template.py \
  --lambdaref-path test/ \
  --paths-md test/MD1/analysis \
  --run-type Debug \
  --plot-type Debug \
```

**Publication mode — two replicas, residue-level convergence, NPZ export:**

```bash
python run_analysis_template.py \
  --lambdaref-path test/ \
  --paths-md test/MD1/analysis test/MD2/analysis \
  --run-type Publication \
  --plot-type Publication \
  --res-ids 75 78 513 \
  --npz-output \
```

**Three-letter labels and custom grid:**

```bash
python run_analysis_template.py \
  --lambdaref-path test/ \
  --paths-md test/MD1/analysis test/MD2/analysis \
  --no-single-letter \
  --plot-rows 15 --plot-cols 4
```

**Homomeric dimer (two chains):**

```bash
python run_analysis_template.py \
  --lambdaref-path test/ \
  --paths-md test/MD1/analysis test/MD2/analysis \
  --chains A B
```

---

## Output Files

All figures are saved in the current working directory. File names are derived from the MD path names.

| File pattern | Generated by | Mode |
|---|---|---|
| `{MD}_histograms.png` | `plot_lambda_hist` | Always |
| `{MD}_timeseries.png` | `plot_protonation_timeseries` | Always |
| `{prefix}_convergence.png` | `plot_protonation_convergence` | Publication |
| `{prefix}_protonfraction.png` | `plot_protonation_fraction` | Publication |
| `Res_{resid}.png` | `single_residue_convergence` | Publication + `--res-ids` |
| `npz_protfrac/*.npz` | statistics export | `--npz-output` |

Where `{prefix}` is automatically derived by joining the parent directory names of `--paths-md` (e.g. `MD1-MD2`).

---

## Module Reference

### analyses.py

**Data loading:**

| Class / Function | Description |
|---|---|
| `XVGData` | Dataclass that loads all `cphmd-coord-*.xvg` files from a directory in parallel using `ProcessPoolExecutor`. Supports indexing by coord ID, slicing, and membership tests. |
| `load_file_for_pool` | Worker function for parallel XVG parsing. Reads time + λ columns, skips comment lines. |

**Statistics:**

| Function | Description |
|---|---|
| `calculate_fractions(array_xvg)` | Classifies frames as protonated (λ < 0.2) or deprotonated (λ > 0.8) and returns both fractions. |
| `get_statistics(coord_id, xvg_data_list)` | Computes mean ± SE of protonation/deprotonation fractions across replicas. |
| `get_protonation_timeseries(coord_id, xvg_data)` | Cumulative protonation fraction as a function of frame index. |

**Histidine-specific** (three λ-coordinates per residue — HSP, HSD, HSE):

| Function | Description |
|---|---|
| `calculate_histidine_fractions(array_xvgs)` | Protonation fraction from the three-state model. |
| `get_histidine_statistics(coord_ids, xvg_data_list)` | Mean ± SE across replicas for histidines. |
| `get_histidine_protonation_timeseries(coord_ids, xvg_data)` | Cumulative protonation time series for histidines. |

**Utilities:**

| Function | Description |
|---|---|
| `parse_lambda_reference(path)` | Parses `lambdareference.dat` into a `{coord_id: lambda_id}` mapping. |
| `resid2coordid(resid, lambda_ref)` | Converts a residue number to its coordinate ID using the lambda reference DataFrame. |

### plot.py

| Function | Returns | Description |
|---|---|---|
| `plot_lambda_hist(...)` | `plt` | Grid of λ-value histograms, one subplot per coordinate. |
| `plot_protonation_timeseries(...)` | `plt` | Grid of cumulative protonation fraction time series. |
| `plot_protonation_convergence(...)` | `plt` | Grid of mean ± SE convergence curves across replicas. |
| `plot_protonation_fraction(...)` | `fig` | Publication-ready bar chart with one axis per residue type. Uses `GridSpec` bin-packing so that all bars have equal physical width and small groups share a row. |
| `single_residue_convergence(...)` | `plt` | Single-panel convergence plot for a specific (non-histidine) residue. |

**Helper functions:**

| Function | Description |
|---|---|
| `_display_resname(resname, single_letter)` | Cleans residue names for labels (strips trailing `T`, converts `HSP` → `HIS`, optionally maps to one-letter codes). |
| `_group_coordids_by_resname(coordids, lambda_ref)` | Groups coordinate IDs by residue type, preserving insertion order. |

---

## Notes on Histidine Handling

Histidine has three protonation micro-states in GROMACS constant-pH, each tracked by its own λ-coordinate:

1. **HSP** (doubly protonated, charged)
2. **HSD** (neutral, proton on Nδ)
3. **HSE** (neutral, proton on Nε)

The toolkit identifies histidines by checking for `resname == "HSPT"` in the lambda reference and always processes three consecutive coordinate IDs as a group. A histidine is counted as **protonated** when the HSP coordinate is > 0.8, and **deprotonated** when either HSD or HSE is > 0.8.

> **Note:** `single_residue_convergence` currently does **not** support histidines. Use the grid-based convergence plot for those.

---

## Multi-Chain (Homomeric) Systems

For homomeric multimers (e.g. a homodimer), the same titratable residue appears in multiple chains with different coordinate IDs. Use `--chains` to enable multi-chain mode:

```bash
python run_analysis_template.py \
  --lambdaref-path project/ \
  --paths-md project/MD1/analysis project/MD2/analysis \
  --chains A B
```

This groups coordinates by chain and builds an internal mapping so that statistics are computed correctly across both chains and replicas. The chain mapping logic may need to be adapted depending on the chain naming convention in your `lambdareference.dat`.

---

## Profiling

To identify performance bottlenecks (useful for large systems or long trajectories), enable cProfile:

```bash
python run_analysis_template.py \
  --lambdaref-path project/ \
  --paths-md project/MD1/analysis \
  --run-type Debug \
  --profile
```

This produces `profile_results.prof` (binary, loadable with `snakeviz` or `pstats`) and `profile_results.txt` (human-readable, sorted by cumulative time).
