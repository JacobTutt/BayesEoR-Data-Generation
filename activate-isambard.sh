#!/usr/bin/env bash
# Source this file from the repository root on Isambard.
module load cray-python/3.11.7 brics/openmpi/4.1.7
export LD_LIBRARY_PATH="${OPENMPI_ROOT}/lib:${LD_LIBRARY_PATH:-}"
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/venv/bin/activate"
