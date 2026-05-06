# Green vs Docking Analysis

`green_vs_docking` scores retrosynthetic routes with green-chemistry metrics, merges the scores with docking rewards, and produces summary tables/plots for condition-optimization experiments.

## Package Layout

Core scoring:
- `run_green_analysis.py`: executable wrapper for the main analysis CLI.
- `main.py`: CLI parsing and output writing for one comparison run.
- `analyzer.py`: loops over experiment folders and seeds.
- `scoring.py`: per-route component scores, green score aggregation, block/protection flags, and summary statistics.
- `atom_economy.py`: route-level atom economy.
- `solvent_score.py`: solvent/agent scoring from `data/solvent_guide.json`.
- `temperature.py`: temperature-range scoring.
- `biomass_utilization.py`: biomass utilisation efficiency (BUE) from block reagents using NameRxn atom mapping.
- `plots.py`: molecule-level green-score/docking merge and KDE plots.

Post-processing:
- `plot_outputs_allexp.py`: aggregate all experiment output folders into one green-vs-docking plot and `resume.csv`.
- `table_summary.py`: reusable logic for the 9KQ3 LaTeX-table summary CSV.
- `build_latex_table_summary.py`: executable wrapper for `table_summary.py`.

Data/configuration:
- `data/blocks.smi`: biomass/block SMILES used for block detection and BUE.
- `data/solvent_guide.json`: solvent/reagent score lookup.
- `environment.yml`: conda environment definition.

## Environment

```bash
conda env create -f /work/liac/tatzber/green_vs_docking/environment.yml
conda activate green-vs-docking
```

## Green Score Definition

For each solved route, the package computes:
- `atom_economy`: route-level heavy-atom mass balance.
- `solvent_score`: minimum solvent/agent score across route steps, using `data/solvent_guide.json` and a default fallback score.
- `temperature_score`: mean score of available reaction-temperature ranges; missing temperatures are ignored.
- `biomass_utilization`: fraction of input biomass block mass that reaches the final product. Routes with no detected block reagent receive `0.0`.

The default green score is an equal-weight average:

```text
green_score = 0.25 atom_economy + 0.25 solvent_score + 0.25 temperature_score + 0.25 biomass_utilization
```

`NaN` components are skipped and the remaining weights are renormalized. Since non-block routes explicitly get `biomass_utilization = 0.0`, the BUE term is included and penalizes routes that do not use blocks.

## Input Contracts

### Experiment Folder

Each experiment path passed to `--experiment` must contain seed folders named `seed0`, `seed1`, etc. Each seed folder must contain:

```text
seedX/
  oracle_history.csv
  syntheseus_results/
    smiles_rxn_tracker.json
```

### `oracle_history.csv`

Required for green/docking analysis:
- `smiles`
- `syntheseus_reward`
- `quickvina2_gpu_raw_values`

A molecule is considered solved by the retrosynthesis model when `syntheseus_reward != 0`.

### `smiles_rxn_tracker.json`

The tracker is expected to map target SMILES to route dictionaries. Route dictionaries should contain `node_*` entries with:
- molecule nodes: `is_mol`, `mol_smiles`, `depth`
- reaction nodes: `is_rxn`, `rxn_smiles`, `depth`, and optionally `agents`, `temperature`, `rxn_class`, `rxn_name`

`rxn_class` and `rxn_name` are used to flag protection/deprotection routes when either field contains the substring `protection`.

### Block and Biomass Inputs

`data/blocks.smi` contains the enforced biomass blocks. Lines can contain one SMILES or comma-separated SMILES. These are canonicalized before matching. BUE requires a NameRxn binary, by default:

```text
/work/liac/tatzber/HazELNut/namerxn
```

## Run Green Analysis

Example:

```bash
python /work/liac/tatzber/green_vs_docking/run_green_analysis.py \
  --experiment baseline=/scratch/sabanza/tango/experiments/conditions/pesticide_basic \
  --experiment pesticide_blocks=/scratch/sabanza/tango/experiments/conditions/pesticide_blocks \
  --output-dir /work/liac/tatzber/green_vs_docking/outputs/pesticide_blocks
```

Optional inputs:

```bash
--solvent-guide /path/to/solvent_guide.json
--block-file /path/to/blocks.smi
--namerxn-binary /path/to/namerxn
```

Outputs written to `--output-dir`:
- `molecule_scores.csv`: one row per scored molecule with component scores, `green_score`, `uses_block`, and `has_protection`.
- `summary_scores.csv`: per-seed rows plus experiment-average rows with mean and standard deviation for each metric.
- `summary_scores_blocks_only.csv`: written when a condition label contains `blocks`; baseline rows are kept unchanged and condition rows are restricted to molecules with `uses_block == True`.
- `green_vs_quickvina.csv`: molecule-level merge of green score and docking value.
- `green_vs_quickvina_kde.png`: KDE plot with docking axis inverted so more negative docking values are on the right.
- `green_vs_quickvina_blocks.csv` / `green_vs_quickvina_blocks_kde.png`: same merge/plot but block-condition rows require `uses_block == True`; baseline remains unchanged.
- `green_vs_quickvina_no_protection.csv` / `green_vs_quickvina_no_protection_kde.png`: same merge/plot but condition rows require `has_protection == False`; baseline remains unchanged.

## Aggregate Existing Experiment Outputs

Use `plot_outputs_allexp.py` when each experiment already has a `green_vs_quickvina.csv` and `summary_scores.csv` output folder.

```bash
python /work/liac/tatzber/green_vs_docking/plot_outputs_allexp.py \
  --root /work/liac/tatzber/green_vs_docking/outputs/9KQ3 \
  --output /work/liac/tatzber/green_vs_docking/outputs/9KQ3/green_vs_docking_summary.png \
  --resume-output /work/liac/tatzber/green_vs_docking/outputs/9KQ3/resume.csv \
  --aggregate-by-seed
```

Useful flags:
- `--aggregate-by-seed`: average each seed first, then compute mean/std across seeds.
- `--prefer-blocks-file`: for folders whose name contains `blocks`, use `green_vs_quickvina_blocks.csv` instead of `green_vs_quickvina.csv`; the default plot filename gets `_blocks` appended if `--output` is not explicitly set.

Default inputs used by this script:
- root: `/work/liac/tatzber/green_vs_docking/outputs/outputs_bae`
- per-experiment comparison file: `green_vs_quickvina.csv`
- per-experiment summaries: `summary_scores.csv` and optional `summary_scores_blocks_only.csv`

Outputs:
- summary scatter plot with docking axis inverted.
- `resume.csv` containing experiment-average metrics from all discovered summary files.

## Build LaTeX-Table Summary CSV

Use this for the 9KQ3 table where molecules are split by experiment and annotated with docking-score, QED, green score, and BMS for solved/enforced molecules.

```bash
python /work/liac/tatzber/green_vs_docking/build_latex_table_summary.py \
  --outputs-root /work/liac/tatzber/green_vs_docking/outputs/9KQ3 \
  --oracle-root /work/liac/vsabanzagil/condition_optimization \
  --output /work/liac/tatzber/green_vs_docking/outputs/9KQ3/latex_table_summary.csv \
  --baseline-experiment pesticide_9KQ3_basic \
  --seeds 0 1 2 3 4
```

The script discovers experiment output folders under `--outputs-root` that contain `green_vs_quickvina.csv`. It matches them to original oracle folders under `--oracle-root` by folder name. Baseline rows are read from the first output folder containing a baseline label, while baseline oracle histories come from `--baseline-experiment`.

Solved definitions for the table:
- all experiments: `syntheseus_reward != 0`
- names containing `blocks`: additionally require `uses_block == True` in `green_vs_quickvina.csv`
- names containing `no_deprotection`: additionally require `has_protection == False` in `green_vs_quickvina.csv`
- if both substrings are present, all conditions must hold

Missing `oracle_history.csv` files are skipped for that experiment. The output records `expected_replicates`, `available_replicates`, `available_seeds`, and `successful_replicates_N` so the table can report exactly which seeds contributed.

Output columns include solved/non-solved mean and standard deviation, pooled solved molecules, docking mean/std, QED mean/std, green-score mean/std, and unique BMS count recomputed from SMILES with RDKit Bemis-Murcko scaffolds.

A Slurm wrapper example is available in `submit_latex_table_summary.sh`.
