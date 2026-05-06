#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=l40s
#SBATCH --gpus=1
#SBATCH --time=2:00:00
#SBATCH --job-name=green_vs_docking

source /work/liac/tatzber/miniconda3/bin/activate tgp

BASELINE="/work/liac/vsabanzagil/condition_optimization/experiment_copy_mateo/pesticide_9KQ3_basic"
CONDITIONS=(
    "/work/liac/vsabanzagil/condition_optimization/experiment_copy_mateo/pesticide_9KQ3_blocks_no_deprotection"
    "/work/liac/vsabanzagil/condition_optimization/experiment_copy_mateo/pesticide_9KQ3_solvent"
    "/work/liac/vsabanzagil/condition_optimization/experiment_copy_mateo/pesticide_9KQ3_solvent_blocks"
    "/work/liac/vsabanzagil/condition_optimization/experiment_copy_mateo/pesticide_9KQ3_solvent_blocks_no_deprotection"
    "/work/liac/vsabanzagil/condition_optimization/experiment_copy_mateo/pesticide_9KQ3_solvent_no_deprotection"
)

for condition in "${CONDITIONS[@]}"; do
    condition_name="$(basename "$condition")"
    python /work/liac/tatzber/green_vs_docking/run_green_analysis.py \
        --experiment baseline="$BASELINE" \
        --experiment "${condition_name}=${condition}" \
        --output-dir "/work/liac/tatzber/green_vs_docking/outputs/9KQ3/${condition_name}"
done
