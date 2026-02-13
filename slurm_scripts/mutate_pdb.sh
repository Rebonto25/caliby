#!/bin/bash   
#SBATCH -J mutate_pdb_fast # Job name
#SBATCH --time=4-00:00:0              # Walltime
#SBATCH --nodes=1                       # number of nodes
#SBATCH --ntasks-per-node=1                      # 1 tasks
#S BATCH --gres=gpu:Quadro_RTX_6000:1      # number of gpus
#S BATCH --gpus-per-task=2               # number of gpus per task
#SBATCH --cpus-per-task=40                # number of cores per task
#SBATCH --mem=20G             # memory/cpu (in MB)
#S BATCH --mem-per-gpu=40000             # memory/gpu (in MB)
#SBATCH --mail-user=rebonto.haque@magd.ox.ac.uk  # set email address
#S BATCH --mail-type=ALL                 # Spam us with everything, caution
#SBATCH --mail-type=begin               # Instead only email when job begins...
#SBATCH --mail-type=end                 # ... and ends
#SBATCH --partition=standard-opig-gpu        # Select a specific partition rather than default
#SBATCH --clusters swan -w nagagpu06.cpu.stats.ox.ac.uk  # Provide a specific node/nodelist rather than the standard no>
#SBATCH --output=/ceph/opig-shared/users/haque/DPhil/caliby/logs/%x/%x_%j.out  # Writes standard output to this file. %j is jobnu>
#SBATCH --error=/ceph/opig-shared/users/haque/DPhil/caliby/logs/%x/%x_error_%j.out   # Writes error messages to this file. %j is >

source ~/.bashrc
conda activate jupyter2
python /ceph/opig-shared/users/haque/DPhil/caliby/util_python_scripts/mutate_pdb.py \
  --parent-pdb /ceph/opig-shared/users/haque/DPhil/inverse_folding/cr9114_H3/H3_CR9114_real_space_refined_003_real_space_refined_005-coot-9_real_space_refined_023.pdb \
  --fasta /ceph/opig-shared/users/haque/DPhil/inverse_folding/cr9114_H3/cr9114_cryo_H3_hc_lib.fasta \
  --out-dir /ceph/opig-shared/users/haque/DPhil/caliby/scoring_analysis/mutant_pdbs_fast \
  --chain H \
  --max-mutations 17 \
  --processes 40 \
  --overwrite