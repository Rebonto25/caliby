#!/bin/bash
# Example script for cleaning PDB files for use with Protpardelle-1c ensemble generation.

source env_setup.sh
python3 caliby/data/preprocessing/atomworks/clean_pdbs.py \
    input_cfg.pdb_dir=/ceph/opig-shared/users/haque/DPhil/caliby/single_structure_pdbs \
    num_workers=50 \
    out_dir=/ceph/opig-shared/users/haque/DPhil/caliby/single_structure_pdbs/cleaned
