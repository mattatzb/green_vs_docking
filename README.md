# Green chemistry vs docking analysis

Reproduce the manuscript's green-chemistry scores, summary tables, and figure
from the archived pesticide experiments. Run the commands below from the
repository root.

## 1. Install the environment

Verified on **Ubuntu 22.04.4 LTS, Linux x86-64**. Analysis dependencies,
including NumPy 1.26.4 and RDKit 2025.3.5, are pinned in `environment.yml`.

```bash
conda env create --file environment.yml
conda activate green-vs-docking
```

Analysis uses saved routes, reaction conditions, and docking scores; it requires
no GPU, model checkpoints, or new docking/retrosynthesis searches.
NameRxn is a separate dependency for route rescoring (step 3).

## 2. Download the archived inputs

Download the required ZIPs from
**[Zenodo: 10.5281/zenodo.23674897](https://doi.org/10.5281/zenodo.23674897)**
(publication pending) and extract them into the repository root:

| Archive | Extracted directories | Purpose |
| --- | --- | --- |
| `green_vs_docking_data.zip` | `data/pesticide_9KQ3_*/` | Saved routes and oracle histories for full rescoring and table generation |
| `green_vs_docking_outputs.zip` | `9KQ3_outputs/pesticide_9KQ3_*/` | Saved per-condition scores, comparisons, summaries, and KDEs for the optional saved-comparison workflow |

```bash
unzip ~/Downloads/green_vs_docking_data.zip -d .
unzip ~/Downloads/green_vs_docking_outputs.zip -d .
```

The data archive contains eight `data/pesticide_9KQ3_*` experiment directories,
including the `pesticide_9KQ3_basic` baseline. Seed inputs follow this layout:

```text
data/pesticide_9KQ3_*/seed0/      # likewise seed1 through seed4
  oracle_history.csv
  syntheseus_results/smiles_rxn_tracker.json
```

`data/blocks.smi`, `data/solvent_guide.json`, and only the three manuscript targets
(`latex_table_summary.csv`, `resume.csv`, and `green_vs_docking_summary_blocks.pdf`)
directly under `9KQ3_outputs/` remain in this repository. The output archive
contains that folder's seven experiment subdirectories. The two root PNGs are
excluded from both Git and the archive. Full rescoring needs only the data archive;
the saved-comparison workflow needs both archives. The two known missing oracle
histories are described below.

## 3. Configure NameRxn for route rescoring

Full route rescoring requires licensed
[NameRxn from NextMove Software](https://www.nextmovesoftware.com/namerxn.html)
for atom mapping. The manuscript used **version 3.7.3**, from
`NameRxn-3.7.3-CentOS-7.8.3-gcc-4.8.5.tar.gz`.
After obtaining and extracting it, configure the executable and check its checksum:

```bash
export NAMERXN_BINARY="/absolute/path/to/HazELNut/namerxn"
sha256sum "$NAMERXN_BINARY"
```

Reference executable SHA-256:
`cd0bc654ab67e633f5537bcacd53eeeb1e302e41a91f771bcf5cc40bd633c43e`.
The saved-comparison workflow below does not require NameRxn.

## 4. Reproduce the full analysis

```bash
python main.py
```

This scores the saved routes for all seven conditions against the baseline,
writes per-condition molecule scores, comparison CSVs, summaries, and KDE plots,
then generates these three manuscript artifacts:

| Generated file | Supplied reference |
| --- | --- |
| `outputs/9KQ3/latex_table_summary.csv` | `9KQ3_outputs/latex_table_summary.csv` |
| `outputs/9KQ3/resume.csv` | `9KQ3_outputs/resume.csv` |
| `outputs/9KQ3/green_vs_docking_summary_blocks.pdf` | `9KQ3_outputs/green_vs_docking_summary_blocks.pdf` |

`python main.py run --output-dir /path/to/fresh-results` selects another
destination. NameRxn can also be selected with `--namerxn-binary`; otherwise
`NAMERXN_BINARY` takes precedence over `namerxn` on `PATH`.
Route scoring and KDE rendering can be slow on the full dataset.

## 5. Check reproduction

Run the small unit and regression tests:

```bash
python -B -m pytest
```

These use synthetic fixtures and require neither the Zenodo downloads nor
NameRxn. They check scoring, statistics options, paths, and pipeline orchestration;
they do not rerun the full manuscript dataset.

After generating the final artifacts, compare both CSVs numerically:

```bash
python - <<'PY'
from pathlib import Path
import pandas as pd

for name in ("latex_table_summary.csv", "resume.csv"):
    pd.testing.assert_frame_equal(
        pd.read_csv(Path("outputs/9KQ3") / name),
        pd.read_csv(Path("9KQ3_outputs") / name),
        check_dtype=False, rtol=1e-10, atol=1e-10,
    )
    print(f"{name}: matches reference")
PY
```

Compare the PDF figure content with the supplied reference; creation metadata
can differ between regenerated PDF files.

## Optional: regenerate from saved comparison CSVs

```bash
python main.py run --from-comparisons 9KQ3_outputs
```

With both ZIPs extracted, this regenerates the same three final artifacts in
`outputs/9KQ3/`, skipping route rescoring and per-condition KDE rendering.
**NameRxn is unnecessary.**
The downloaded oracle histories are still required for the LaTeX-table summary.
The reference inputs remain unchanged with the default destinations.

## Reproducibility limits and implementation notes

- **Missing histories:** interrupted SATURN runs left no `oracle_history.csv`
  for `pesticide_9KQ3_blocks/seed3` and
  `pesticide_9KQ3_blocks_no_deprotection/seed4`. These histories cannot be recovered.
  Seeds lacking required inputs are skipped; the table records
  `expected_replicates`, `available_replicates`, `available_seeds`, and
  `successful_replicates_N`.
- **Green score:** the equal-weight average of atom economy, solvent score,
  temperature score, and biomass utilization uses the supplied block list and
  solvent guide. Missing components are omitted and remaining weights renormalized.
- **Manuscript statistics:** the plot and `resume.csv` pool molecule rows across
  seeds, report sample standard deviations, and select
  `green_vs_quickvina_blocks.csv` for block conditions. The table reports
  solved/non-solved counts across seeds and molecular metrics over the pooled
  solved set. Solved molecules require nonzero `syntheseus_reward` and an available
  green score; block conditions additionally require block use, and
  `no_deprotection` conditions require no protection/deprotection flag.
- **Alternative statistics:** `--aggregate-by-seed` averages seed means for the
  plot and resume; `--no-prefer-blocks-file` selects unfiltered comparisons.
  `--resume-from-summaries` exports historical summary-file averages to the resume
  while retaining pooled plot statistics, and cannot be combined with
  `--aggregate-by-seed`. These options do not change the table's definitions.
- **Architecture:** `main.py` dispatches the `run`, `score`, `plot`, and `table`
  commands. `src/green_vs_docking/analysis/` computes metrics and tables,
  `plots/` renders figures, and `cli/` provides command-line interfaces.
  See `python main.py <command> --help` for stage options, including optional
  joint KDE plots. `pyproject.toml` configures packaging and the public test suite;
  running `main.py` does not require installing the package.
