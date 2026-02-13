#!/bin/bash
set -euo pipefail

# Run score_ensemble.py for each mutant PDB by reusing a single ensemble directory.
# The script renames the ensemble directory to the mutant ID, replaces the primary
# PDB file to match that ID, then runs scoring.
#
# Usage:
#   examples/scripts/score_mutant_ensembles.sh \
#     --mutants /path/to/mutant_pdbs \
#     --ensemble-dir /path/to/generic_run_name/ORIGINAL_PDB_ID \
#     --ckpt model_params/caliby/caliby.ckpt \
#     --out-root examples/outputs/score_ensemble_mutants
#
# Notes:
# - Mutant PDB filenames are assumed to be <PDB_ID>.pdb (stem used as ID).
# - Conformer files inside the ensemble dir can keep their existing names.
# - The primary conformer must be named <PDB_ID>.pdb, and the directory must be named <PDB_ID>.

source env_setup.sh

MUTANTS_DIR=""
ENSEMBLE_DIR=""
CKPT_PATH=""
OUT_ROOT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mutants)
      MUTANTS_DIR="$2"; shift 2 ;;
    --ensemble-dir)
      ENSEMBLE_DIR="$2"; shift 2 ;;
    --ckpt)
      CKPT_PATH="$2"; shift 2 ;;
    --out-root)
      OUT_ROOT="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1"; exit 1 ;;
  esac
done

if [[ -z "$MUTANTS_DIR" || -z "$ENSEMBLE_DIR" || -z "$CKPT_PATH" || -z "$OUT_ROOT" ]]; then
  echo "Missing required arguments."
  exit 1
fi

if [[ ! -d "$MUTANTS_DIR" ]]; then
  echo "Mutants directory not found: $MUTANTS_DIR"; exit 1
fi

if [[ ! -d "$ENSEMBLE_DIR" ]]; then
  echo "Ensemble directory not found: $ENSEMBLE_DIR"; exit 1
fi

PARENT_DIR="$(dirname "$ENSEMBLE_DIR")"
CURRENT_DIR_NAME="$(basename "$ENSEMBLE_DIR")"

shopt -s nullglob
for pdb in "$MUTANTS_DIR"/*.pdb; do
  mutant_id="$(basename "$pdb" .pdb)"

  # Rename ensemble dir to match mutant ID (required by scoring code)
  if [[ "$CURRENT_DIR_NAME" != "$mutant_id" ]]; then
    mv "$PARENT_DIR/$CURRENT_DIR_NAME" "$PARENT_DIR/$mutant_id"
    CURRENT_DIR_NAME="$mutant_id"
  fi

  # Replace primary conformer with the mutant PDB
  rm -f "$PARENT_DIR/$CURRENT_DIR_NAME/$CURRENT_DIR_NAME.pdb"
  cp "$pdb" "$PARENT_DIR/$CURRENT_DIR_NAME/$CURRENT_DIR_NAME.pdb"

  # Run scoring for this mutant
  python3 caliby/eval/sampling/score_ensemble.py \
    ckpt_path="$CKPT_PATH" \
    input_cfg.conformer_dir="$PARENT_DIR" \
    input_cfg.pdb_name_list="[$mutant_id]" \
    out_dir="$OUT_ROOT/$mutant_id"

done
