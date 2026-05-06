#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=8
#SBATCH --partition=l40s
#SBATCH --gpus=1
#SBATCH --job-name=9KQ3_table 
#SBATCH --time=02:00:00

source /work/liac/tatzber/miniconda3/etc/profile.d/conda.sh
conda activate tgp

python /work/liac/tatzber/green_vs_docking/build_latex_table_summary.py \
  --outputs-root /work/liac/tatzber/green_vs_docking/outputs/9KQ3 \
  --oracle-root /work/liac/vsabanzagil/condition_optimization \
  --output /work/liac/tatzber/green_vs_docking/outputs/9KQ3/latex_table_summary.csv
