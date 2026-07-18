#!/usr/bin/env bash
#SBATCH --job-name=beor-generate-skies
#SBATCH --account=u6my
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --mem=120G
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --output=logs/prepare-%j.out

set -euo pipefail
repo_dir="${SLURM_SUBMIT_DIR}"
mkdir -p "${repo_dir}/logs"
module load cray-python/3.11.7
source "${repo_dir}/venv/bin/activate"
python "${repo_dir}/generate.py" "$@" prepare
