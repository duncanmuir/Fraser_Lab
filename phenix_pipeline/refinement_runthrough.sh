#!/bin/bash
# set up a new refinement
# Stephanie Wankowicz 4/5/2019
# updated: 11-01-2019
# modified by Duncan Muir 11/8/2019

# source Phenix
source /Applications/phenix-1.16-3549/phenix_env.sh
export PHENIX_OVERWRITE_ALL=true

PDB=$1

NSLOTS=$2
echo $PDB
echo "$PWD"

# if you want to run amber
# phenix.AmberPrep ${PDB}.pdb

# get ligand name and list of ligands for harmonic restraints
~/opt/anaconda3/bin/python ~/Fraser_Lab/phenix/PDB_ligand_parser.py  ~/Fraser_Lab/phenix_pipeline/ligands_to_remove.csv
lig_name=$(cat "ligand_name.txt")
echo "$lig_name"

echo '________________________________________________________Starting Phenix elbow________________________________________________________'
#phenix.elbow $PDB.pdb --residue $lig_name --final_geometry

echo '________________________________________________________Starting Phenix cif as mtz________________________________________________________'
phenix.cif_as_mtz "$PDB"-sf.cif --extend_flags --merge

echo '________________________________________________________Starting Phenix Ready Set________________________________________________________'
phenix.ready_set pdb_file_name="$PDB".pdb #cif_file_name=elbow.${lig_name}.${PDB}_pdb.001.cif >> readyset_output.txt

echo '________________________________________________________Checking on FOBS________________________________________________________'
if grep -F _refln.F_meas_au "$PDB"-sf.cif; then
        echo 'FOBS'
else
        echo 'IOBS'
fi
rm "${PDB}".updated_refine_*
if [[ -e "4phenix_${PDB}.ligands.cif" ]]; then
echo '________________________________________________________Running refinement with ligand.________________________________________________________'
    if grep -F _refln.F_meas_au "$PDB"-sf.cif; then
        phenix.refine "$PDB".updated.pdb "$PDB"-sf.mtz "${PDB}".ligands.cif refinement.input.xray_data.r_free_flags.label=R-free-flags ~/Fraser_Lab/phenix_pipeline/finalize.params refinement.input.xray_data.labels="FOBS,SIGFOBS" # nproc=$NSLOTS use_amber=True amber.topology_file_name=4amber_${PDB}.prmtop amber.coordinate_file_name=4amber_${PDB}.rst7 amber.order_file_name=4amber_${PDB}.order
    else
        echo 'IOBS'   
         phenix.refine "$PDB".updated.pdb "$PDB"-sf.mtz 4phenix_"${PDB}".ligands.cif refinement.input.xray_data.r_free_flags.label=R-free-flags  refinement.input.xray_data.labels="IOBS,SIGIOBS" # nproc=$NSLOTS use_amber=True amber.topology_file_name=4amber_${PDB}.prmtop amber.coordinate_file_name=4amber_${PDB}.rst7 amber.order_file_name=4amber_${PDB}.order 4phenix_$PDB.updated.pdb
    fi
    
else
    echo '________________________________________________________Running refinement without ligand.________________________________________________________'
    if grep -F _refln.F_meas_au "$PDB"-sf.cif; then
        phenix.refine "$PDB".updated.pdb "$PDB"-sf.mtz ~/Fraser_Lab/phenix_pipeline/finalize.params refinement.input.xray_data.labels="FOBS,SIGFOBS" refinement.input.xray_data.r_free_flags.label=R-free-flags # nproc=$NSLOTS use_amber=True topology_file_name=${PDB}.prmtop coordinate_file_name=4amber_${PDB}.rst7
    else
        phenix.refine "$PDB".updated.pdb "$PDB"-sf.mtz ~/Fraser_Lab/phenix_pipeline/finalize.params refinement.input.xray_data.r_free_flags.label=R-free-flags refinement.input.xray_data.labels="IOBS,SIGIOBS" # nproc=$NSLOTS use_amber=True topology_file_name=${PDB}.prmtop coordinate_file_name=${PDB}.rst7
   fi
 fi

echo '________________________________________________________Starting Model versus Data________________________________________________________'
phenix.model_vs_data "${PDB}".updated_refine_001.pdb $PDB-sf.cif > "${PDB}"_model_v_data.txt

echo '________________________________________________________Starting extract python script________________________________________________________'
  ~/opt/anaconda3/bin/python ~/Fraser_Lab/phenix/parse_refine_log.py  -log_file "${PDB}".updated_refine_001.log -PDB $PDB
  rm lig_RMSZ_updated.txt
  rm lig_RMSZ_pre_refine.txt
  
  echo '________________________________________________________Validating the Ligand from Original PDB________________________________________________________'
  if [[ -e "${PDB}.ligands.cif" ]]; then  #only running this if we have a ligand    
    
    echo "$lig_name"
    #mmtbx.validate_ligands ${PDB}.pdb ${PDB}-sf.mtz ligand_code=$lig_name #"${PDB}_ligand" #prints out ADPs and occs + additional information
    
    echo '________________________________________________________Phenix PDB Interpretation_______________________________________________________'
    phenix.pdb_interpretation "${PDB}".updated.pdb "${PDB}".ligands.cif write_geo=True
    
    echo '________________________________________________________Elbow Refine_Geo_Display_______________________________________________________'
    elbow.refine_geo_display "${PDB}".updated.pdb.geo "$lig_name" >> lig_RMSZ_pre_refine.txt #prints out deviations including ligand specific RMSD and RMSz values
  fi


echo '________________________________________________________Validating the Ligand after Initial Refinement________________________________________________________'
  if [[ -e "${PDB}.ligands.cif" ]]; then  #only running this if we have a ligand
    #mmtbx.validate_ligands ${PDB}.updated_refine_001.pdb ${PDB}-sf.mtz ligand_code=$lig_name #"${PDB}_ligand" #prints out ADPs and occs + additional information
    
    echo '________________________________________________________Phenix PDB Interpretation_______________________________________________________'
    phenix.pdb_interpretation "${PDB}".updated_refine_001.pdb "${PDB}".ligands.cif write_geo=True
    
    echo '________________________________________________________Elbow Refine_Geo_Display_______________________________________________________'
    elbow.refine_geo_display "${PDB}".updated_refine_001.pdb.geo "$lig_name" >> lig_RMSZ_updated.txt
    #prints out deviations including ligand specific RMSD and RMSz values
    #calculate the energy difference of the ligand in the model and it relaxed RM1/AM1, but can be linked to 3rd party packages
  fi

  echo '________________________________________________________Validating Ligand Output_______________________________________________________'
  ~/opt/anaconda3/bin/python ~/Fraser_Lab/phenix/ligand_geo_parser.py -pre_refine=lig_RMSZ_pre_refine.txt -post_refine=lig_RMSZ_updated.txt -PDB=$PDB


#  echo '________________________________________________________Begin Ensemble Refinement_______________________________________________________'
#  qsub /wynton/home/fraserlab/swankowicz/190419_Phenix_ensemble/grid_search_ens_refine.sh $PDB $3

