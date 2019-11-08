#!/bin/bash
# Stephanie Wankowicz
# 04/25/2019
# Modified by Duncan Muir, 11/08/2019

# this must be done before you submit to SGE since SGE cannot connect to the internet!

#________________________________________________INPUTS________________________________________________
base_folder='/Users/dmuir/kinome_ensemble_refinement/initial_egfr_refinements' #base folder

pdb_filelist=~/kinome_ensemble_refinement/initial_egfr_refinements/pdb_list.txt #list of pdb files
while read -r line; do
  PDB=$line
  cd $base_folder
  if [ -d "/$PDB" ]; then
    echo "Folder exists."
  else
    mkdir $PDB
  fi
  #mkdir $PDB
  cd $PDB
#  phenix.fetch_pdb $PDB
#  phenix.fetch_pdb $PDB --mtz
#  phenix.fetch_pdb $PDB -x
  wget https://files.rcsb.org/download/${PDB}.pdb
  wget https://files.rcsb.org/download/${PDB}-sf.cif
  wget http://edmaps.rcsb.org/coefficients/${PDB}.mtz
done < $pdb_filelist





