#!/bin/bash
# Example script for generating ensembles with Protpardelle-1c.

source env_setup.sh
python3 caliby/eval/sampling/generate_ensembles.py \
    model_params_path=model_params \
    input_cfg.pdb_dir=/ceph/opig-shared/users/haque/DPhil/caliby/single_structure_pdbs \
    num_samples_per_pdb=32 \
    out_dir=/ceph/opig-shared/users/haque/DPhil/caliby/ensemble_pdbs
