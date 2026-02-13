#!/bin/bash   
#SBATCH -J generate_ensembles_H3 # Job name
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
#SBATCH --output=/ceph/opig-shared/users/haque/DPhil/logs/%x/%x_%j.out  # Writes standard output to this file. %j is jobnu>
#SBATCH --error=/ceph/opig-shared/users/haque/DPhil/logs/%x/%x_error_%j.out   # Writes error messages to this file. %j is >

cd /ceph/opig-shared/users/haque/DPhil/caliby
source env_setup.sh
python3 caliby/eval/sampling/generate_ensembles.py \
    model_params_path=model_params \
    input_cfg.pdb_dir=/ceph/opig-shared/users/haque/DPhil/caliby/cleaned_single_structure_pdb/cryo_H3 \
    num_samples_per_pdb=32 \
    out_dir=/ceph/opig-shared/users/haque/DPhil/caliby/ensemble_pdbs/cleaned_ensemble_pdbs/cryo_H3
