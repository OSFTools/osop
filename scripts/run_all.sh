#!/usr/bin/env bash
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG_FILE="${REPO_ROOT}/osop_config.yml"

# Legacy compatibility: map -t to the new --test-mode flag.
if [[ "${1:-}" == "-t" ]]; then
    shift
    exec python "${SCRIPT_DIR}/run_all.py" --config "${CONFIG_FILE}" --test-mode "$@"
fi

exec python "${SCRIPT_DIR}/run_all.py" --config "${CONFIG_FILE}" "$@"
