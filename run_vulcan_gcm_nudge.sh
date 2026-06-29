#!/usr/bin/env bash

set -euo pipefail

VULCAN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_DIR="$(cd "${VULCAN_DIR}/.." && pwd)"

CONDA_BASE="${CONDA_BASE:-/Users/haolinli/miniconda3}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-Vulcan}"
CONDA_SH="${CONDA_BASE}/etc/profile.d/conda.sh"

GCM_INPUT="${GCM_INPUT:-${WORKSPACE_DIR}/GCM_data/day0070h00.atmos_ins.nc}"
VULCAN_CFG="${VULCAN_CFG:-${VULCAN_DIR}/vulcan_cfg.py}"
GCM_OUTPUT_DIR="${GCM_OUTPUT_DIR:-${VULCAN_DIR}/atm/gcm_columns}"
BACKGROUND_MMW="${BACKGROUND_MMW:-44.011}"
GCM_STEM="${GCM_STEM:-$(basename "${GCM_INPUT%.nc}")}"
HYBRID_COEFF_FILE="${HYBRID_COEFF_FILE:-${WORKSPACE_DIR}/GCM_data/${GCM_STEM}_hybrid_coeffs.nc}"
HYBRID_PS0_PA="${HYBRID_PS0_PA:-1000000.0}"
HYBRID_PTOP_PA="${HYBRID_PTOP_PA:-100.0}"
HYBRID_PINT_FACTOR="${HYBRID_PINT_FACTOR:-0.4}"
HYBRID_REFERENCE_PRESSURE_NC="${HYBRID_REFERENCE_PRESSURE_NC:-}"

if [[ ! -f "${CONDA_SH}" ]]; then
  echo "Cannot find conda activation script: ${CONDA_SH}" >&2
  exit 1
fi

if [[ ! -f "${GCM_INPUT}" ]]; then
  echo "Cannot find GCM input file: ${GCM_INPUT}" >&2
  exit 1
fi

if [[ ! -f "${VULCAN_CFG}" ]]; then
  echo "Cannot find VULCAN config: ${VULCAN_CFG}" >&2
  exit 1
fi

source "${CONDA_SH}"
conda activate "${CONDA_ENV_NAME}"

export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mpl-vulcan}"

mkdir -p "${VULCAN_DIR}/output" "${VULCAN_DIR}/plot" "${VULCAN_DIR}/plot/movie" "${GCM_OUTPUT_DIR}"

echo "==> Building hybrid coefficient sidecar"
HYBRID_CMD=(
  python "${VULCAN_DIR}/tools/make_uniform_hybrid_coeffs.py"
  "${GCM_INPUT}"
  --output-path "${HYBRID_COEFF_FILE}"
  --ps0-pa "${HYBRID_PS0_PA}"
  --ptop-pa "${HYBRID_PTOP_PA}"
  --pint-factor "${HYBRID_PINT_FACTOR}"
)
if [[ -n "${HYBRID_REFERENCE_PRESSURE_NC}" ]]; then
  HYBRID_CMD+=(--reference-pressure-nc "${HYBRID_REFERENCE_PRESSURE_NC}")
fi
"${HYBRID_CMD[@]}"

echo "==> Rebuilding GCM nudging driver"
python "${VULCAN_DIR}/tools/extract_gcm_column_for_vulcan.py" \
  "${GCM_INPUT}" \
  --hybrid-coeff-nc "${HYBRID_COEFF_FILE}" \
  --background-mean-mol-weight "${BACKGROUND_MMW}" \
  --vulcan-cfg "${VULCAN_CFG}" \
  --output-dir "${GCM_OUTPUT_DIR}"

echo "==> Validating chemical network"
cd "${VULCAN_DIR}"
python -u make_chem_funs.py

echo "==> Running VULCAN"
python -u vulcan.py -n

echo "==> Finished"
echo "Output: ${VULCAN_DIR}/output/HD189-gcm-nudge.vul"
echo "Plots:  ${VULCAN_DIR}/plot/HD189-gcm-nudge_time_series.png"
echo "Plots:  ${VULCAN_DIR}/plot/HD189-gcm-nudge_profiles.png"
