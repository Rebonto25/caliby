#!/bin/bash
#SBATCH -J rebuild_sidechains
#SBATCH --time=4-00:00:0              # Walltime
#SBATCH --nodes=1                       # number of nodes
#SBATCH --ntasks-per-node=1                      # 1 tasks
#S BATCH --gres=gpu:Quadro_RTX_6000:1      # number of gpus
#S BATCH --gpus-per-task=2               # number of gpus per task
#SBATCH --cpus-per-task=45                # number of cores per task
#SBATCH --mem=20G             # memory/cpu (in MB)
#S BATCH --mem-per-gpu=40000             # memory/gpu (in MB)
#SBATCH --mail-user=rebonto.haque@magd.ox.ac.uk  # set email address
#S BATCH --mail-type=ALL                 # Spam us with everything, caution
#SBATCH --mail-type=begin               # Instead only email when job begins...
#SBATCH --mail-type=end                 # ... and ends
#SBATCH --partition=standard-opig-gpu        # Select a specific partition rather than default
#SBATCH --clusters=swan
#SBATCH --nodelist=nagagpu06.cpu.stats.ox.ac.uk
#SBATCH --output=/ceph/opig-shared/users/haque/DPhil/caliby/logs/%x/%x_%A_%a.out
#SBATCH --error=/ceph/opig-shared/users/haque/DPhil/caliby/logs/%x/%x_error_%A_%a.out
set -euo pipefail
cd /ceph/opig-shared/users/haque/DPhil/caliby
source env_setup.sh
python util_python_scripts/rebuild_sidechains_simple.py --input-dir /ceph/opig-shared/users/haque/DPhil/caliby/scoring_analysis/mutant_pdbs_fast --output-dir /ceph/opig-shared/users/haque/DPhil/caliby/scoring_analysis/mutant_pdbs_fullatom_2 --workers 45