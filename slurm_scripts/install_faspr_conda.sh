#!/bin/bash   
#SBATCH -J install_faspr # Job name
#SBATCH --time=4-00:00:0              # Walltime
#SBATCH --nodes=1                       # number of nodes
#SBATCH --ntasks-per-node=1                      # 1 tasks
#SBATCH --gres=gpu:Ada_RTX_6000:1      # number of gpus
#S BATCH --gpus-per-task=2               # number of gpus per task
#SBATCH --cpus-per-task=1                # number of cores per task
#SBATCH --mem=20G             # memory/cpu (in MB)
#S BATCH --mem-per-gpu=40000             # memory/gpu (in MB)
#SBATCH --mail-user=rebonto.haque@magd.ox.ac.uk  # set email address
#S BATCH --mail-type=ALL                 # Spam us with everything, caution
#SBATCH --mail-type=begin               # Instead only email when job begins...
#SBATCH --mail-type=end                 # ... and ends
#SBATCH --partition=standard-opig-gpu        # Select a specific partition rather than default
#SBATCH --clusters swan -w nagagpu09.cpu.stats.ox.ac.uk  # Provide a specific node/nodelist rather than the standard no>
#SBATCH --output=/ceph/opig-shared/users/haque/DPhil/caliby/logs/%x/%x_%j.out  # Writes standard output to this file. %j is jobnu>
#SBATCH --error=/ceph/opig-shared/users/haque/DPhil/caliby/logs/%x/%x_error_%j.out   # Writes error messages to this file. %j is >

set -euo pipefail

# Vanilla FASPR installer (matches upstream GitHub instructions).
#
# Usage:
#   bash slurm_scripts/install_faspr_conda.sh
#   bash slurm_scripts/install_faspr_conda.sh --src-dir "$HOME/src/FASPR"
#   bash slurm_scripts/install_faspr_conda.sh --install-dir "$HOME/.local/opt/faspr"
#   bash slurm_scripts/install_faspr_conda.sh --repo https://github.com/tommyhuangthu/FASPR.git

SRC_DIR="/ceph/opig-shared/users/haque/DPhil/FASPR"
INSTALL_DIR="$HOME/.local/opt/faspr"
REPO_URL="https://github.com/tommyhuangthu/FASPR.git"
JOBS="$(nproc)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --src-dir)
      SRC_DIR="$2"
      shift 2
      ;;
    --install-dir)
      INSTALL_DIR="$2"
      shift 2
      ;;
    --prefix)
      INSTALL_DIR="$2"
      shift 2
      ;;
    --repo)
      REPO_URL="$2"
      shift 2
      ;;
    --jobs)
      JOBS="$2"
      shift 2
      ;;
    -h|--help)
      sed -n '1,45p' "$0"
      exit 0
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

for cmd in git g++; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: required command not found: $cmd"
    exit 1
  fi
done

mkdir -p "$SRC_DIR" "$INSTALL_DIR" "$HOME/.local/bin"

if [[ ! -d "$SRC_DIR/.git" ]]; then
  echo "[1/5] Cloning FASPR..."
  git clone "$REPO_URL" "$SRC_DIR"
else
  echo "[1/5] Updating FASPR..."
  git -C "$SRC_DIR" pull --ff-only
fi

if [[ ! -d "$SRC_DIR/src" ]]; then
  echo "ERROR: expected source directory not found: $SRC_DIR/src"
  exit 1
fi

echo "[2/5] Building with g++ (upstream method)..."
g++ -O3 --fast-math -o "$SRC_DIR/FASPR" "$SRC_DIR"/src/*.cpp

if [[ ! -f "$SRC_DIR/dun2010bbdep.bin" ]]; then
  echo "ERROR: dun2010bbdep.bin not found in source tree."
  echo "FASPR requires dun2010bbdep.bin to be in the same directory as the executable."
  exit 1
fi

echo "[3/5] Installing executable and rotamer library..."
cp -f "$SRC_DIR/FASPR" "$INSTALL_DIR/FASPR"
cp -f "$SRC_DIR/dun2010bbdep.bin" "$INSTALL_DIR/dun2010bbdep.bin"

echo "[4/5] Creating user-level launcher..."
cat > "$HOME/.local/bin/faspr" <<EOF
#!/usr/bin/env bash
exec "$INSTALL_DIR/FASPR" "\$@"
EOF
chmod +x "$HOME/.local/bin/faspr"

echo "[5/5] Final checks..."
if [[ ! -x "$INSTALL_DIR/FASPR" ]]; then
  echo "ERROR: installed executable not found: $INSTALL_DIR/FASPR"
  exit 1
fi

if [[ ! -f "$INSTALL_DIR/dun2010bbdep.bin" ]]; then
  echo "ERROR: rotamer library not found: $INSTALL_DIR/dun2010bbdep.bin"
  exit 1
fi

if [[ -x "$HOME/.local/bin/faspr" ]]; then
  BIN_PATH="$HOME/.local/bin/faspr"
else
  echo "ERROR: FASPR executable not found after install."
  exit 1
fi

echo ""
echo "FASPR is ready for fast side-chain packing."
echo "Binary: $BIN_PATH"
echo "Install dir (contains both files): $INSTALL_DIR"
echo "  - $INSTALL_DIR/FASPR"
echo "  - $INSTALL_DIR/dun2010bbdep.bin"
echo "Version/help: $INSTALL_DIR/FASPR -h"
