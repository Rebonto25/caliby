#!/bin/bash   
#SBATCH -J score_mutants_ensemble_test # Job name
#SBATCH --time=4-00:00:0              # Walltime
#SBATCH --nodes=1                       # number of nodes
#SBATCH --ntasks-per-node=1                      # 1 tasks
#SBATCH --gres=gpu:Ada_RTX_6000:1      # number of gpus
#S BATCH --gpus-per-task=2               # number of gpus per task
#SBATCH --cpus-per-task=1                # number of cores per task
#SBATCH --mem=50G             # memory/cpu (in MB)
#S BATCH --mem-per-gpu=40000             # memory/gpu (in MB)
#SBATCH --mail-user=rebonto.haque@magd.ox.ac.uk  # set email address
#S BATCH --mail-type=ALL                 # Spam us with everything, caution
#SBATCH --mail-type=begin               # Instead only email when job begins...
#SBATCH --mail-type=end                 # ... and ends
#SBATCH --partition=standard-opig-gpu        # Select a specific partition rather than default
#SBATCH --clusters swan -w nagagpu09.cpu.stats.ox.ac.uk  # Provide a specific node/nodelist rather than the standard no>
#SBATCH --output=/ceph/opig-shared/users/haque/DPhil/caliby/logs/%x/%x_%j.out  # Writes standard output to this file. %j is jobnu>
#SBATCH --error=/ceph/opig-shared/users/haque/DPhil/caliby/logs/%x/%x_error_%j.out   # Writes error messages to this file. %j is >

cd /ceph/opig-shared/users/haque/DPhil/caliby
source env_setup.sh

# OOM mitigation: reduce fragmentation and silence unrelated deprecation spam.
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"
export PYTHONWARNINGS="ignore::DeprecationWarning"

# Main memory lever in ensemble scoring is number of conformers loaded at once.
MAX_NUM_CONFORMERS="2"
INCLUDE_PRIMARY_CONFORMER="true"

python3 caliby/eval/sampling/score_ensemble.py \
    ckpt_path=model_params/caliby/caliby.ckpt \
    input_cfg.conformer_dir=/ceph/opig-shared/users/haque/DPhil/caliby/ensemble_pdbs/cleaned_ensemble_pdbs/cryo_H3/cc95-epoch3490-sampling_partial_diffusion-ss1.0-schurn0-ccstart0.0-dx0.0-dy0.0-dz0.0-rewind150 \
    out_dir=examples/outputs/score_ensemble \
    max_num_conformers=${MAX_NUM_CONFORMERS} \

