#!/usr/bin/env bash
#SBATCH --job-name=bayeseor-vis
#SBATCH --account=brics.u6my
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=18
#SBATCH --mem=120G
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=logs/simulate-%j.out

set -euo pipefail

repository="${SLURM_SUBMIT_DIR}"
mkdir -p "${repository}/logs"
module load cray-python/3.11.7 brics/openmpi/4.1.7
export LD_LIBRARY_PATH="${OPENMPI_ROOT}/lib:${LD_LIBRARY_PATH:-}"
source "${repository}/venv/bin/activate"

srun python "${repository}/simulate.py" "$@"
