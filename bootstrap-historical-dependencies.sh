#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dependency_dir="${repo_dir}/external/BayesEoR"
commit="86d5f9239b46610d5226d699a1a0c2e2602c5ea5"

if [[ ! -d "${dependency_dir}/.git" ]]; then
  git clone https://github.com/PSims/BayesEoR.git "${dependency_dir}"
fi

git -C "${dependency_dir}" fetch origin "${commit}"
git -C "${dependency_dir}" checkout --detach "${commit}"
echo "Historical BayesEoR preprocessor pinned at ${commit}"
