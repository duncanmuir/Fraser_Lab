#!/usr/bin/env python
import subprocess as sp
import argparse
import wget
import os
import sys
import shutil
from biopandas.pdb import PandasPdb
import pandas as pd


def yes_or_no(question):
    while "the answer is invalid":
        reply = str(input(question + " (y/n): ".lower().strip()))
        if reply[:1] == 'y':
            return True
        if reply[:1] == 'n':
            return False


def download_input_files(pdb_id, base_folder):
    """
    Fetches refinement input files from RCSB PDB
    :param pdb_id:
    :return:
    """

    output_path = base_folder + "/" + pdb_id

    if os.path.exists(output_path):
        if yes_or_no("Directory already exists, would you like to remove?"):
            shutil.rmtree(output_path)
        else:
            sys.exit(0)

    os.mkdir(output_path)

    wget.download("https://files.rcsb.org/download/{}.pdb".format(pdb_id),
                  "{}/{}.pdb".format(output_path, pdb_id))
    wget.download("https://files.rcsb.org/download/{}-sf.cif".format(pdb_id),
                  "{}/{}-sf.cif".format(output_path, pdb_id))
    wget.download("http://edmaps.rcsb.org/coefficients/{}.mtz".format(pdb_id),
                  "{}/{}.mtz".format(output_path, pdb_id))

# if you want to run amber
# phenix.AmberPrep ${PDB}.pdb


def parse_ligand_from_pdb(pdb_id, base_folder):
    """
    Identifies drug-like ligands from PDB input file.
    :param pdb_id:
    :param base_folder:
    :return:
    """

    # Read PDB file into PandasPDB df
    ppdb = PandasPdb()
    ppdb.read_pdb("{}/{}/{}.pdb".format(base_folder,pdb_id,pdb_id))

    # Subset df to hetatms
    hetatm_df = ppdb.df['HETATM']

    # Read in ligands
    lig_to_remove_df = pd.read_csv("~/Fraser_Lab/phenix_pipeline/ligands_to_remove.csv")
    lig_to_remove_df.columns = ["name", "unknown"]

    # Get list of unique residue names
    residue_names = list(set(hetatm_df['residue_name']))

    lig_name = ''
    for res_name in residue_names:
        # all_ligands.append(i)
        if res_name not in set(lig_to_remove_df['name']):
            print("###################################")
            lig_name = res_name
            print(lig_name)

    return lig_name


def run_refinement(pdb_id):
    if os.path.exists(f'{pdb_id}.ligands.cif'):

        print("\nRunning refinement with ligand...\n")
        if sp.check_output(f'grep -F _refln.F_meas_au {pdb_id}-sf.cif', shell=True).decode("utf-8") != "":
            sp.call(f'phenix.refine {pdb_id}.updated.pdb {pdb_id}-sf.mtz {pdb_id}.ligands.cif ' 
                    f'refinement.input.xray_data.r_free_flags.label=R-free-flags '
                    f'~/Fraser_Lab/phenix_pipeline/finalize.params refinement.input.xray_data.labels="FOBS,SIGFOBS"',
                    shell=True)
            # nproc=$NSLOTS use_amber=True amber.topology_file_name=4amber_${PDB}.prmtop amber.coordinate_file_name=4amber_${PDB}.rst7 amber.order_file_name=4amber_${PDB}.order
        else:
            'IOBS'
            sp.call(f'phenix.refine {pdb_id}.updated.pdb {pdb_id}-sf.mtz {pdb_id}.ligands.cif'
                    f'refinement.input.xray_data.r_free_flags.label=R-free-flags'
                    f'~/Fraser_Lab/phenix_pipeline/finalize.params refinement.input.xray_data.labels="IOBS,SIGIOBS"',
                    shell=True)
                    # nproc=$NSLOTS use_amber=True amber.topology_file_name=4amber_${PDB}.prmtop amber.coordinate_file_name=4amber_${PDB}.rst7 amber.order_file_name=4amber_${PDB}.order 4phenix_$PDB.updated.pdb

    else:
        print("\nRunning refinement without ligand...\n")
        if sp.check_output(f'grep -F _refln.F_meas_au {pdb_id}-sf.cif', shell=True).decode("utf-8") != "":
            sp.call(f'phenix.refine {pdb_id}.updated.pdb {pdb_id}-sf.mtz '
                    f'refinement.input.xray_data.r_free_flags.label=R-free-flags '
                    f'~/Fraser_Lab/phenix_pipeline/finalize.params refinement.input.xray_data.labels="FOBS,SIGFOBS"',
                    shell=True)
            # nproc=$NSLOTS use_amber=True amber.topology_file_name=4amber_${PDB}.prmtop amber.coordinate_file_name=4amber_${PDB}.rst7 amber.order_file_name=4amber_${PDB}.order
        else:
            'IOBS'
            sp.call(f'phenix.refine {pdb_id}.updated.pdb {pdb_id}-sf.mtz '
                    f'refinement.input.xray_data.r_free_flags.label=R-free-flags'
                    f'~/Fraser_Lab/phenix_pipeline/finalize.params refinement.input.xray_data.labels="IOBS,SIGIOBS"',
                    shell=True)
            # nproc=$NSLOTS use_amber=True amber.topology_file_name=4amber_${PDB}.prmtop amber.coordinate_file_name=4amber_${PDB}.rst7 amber.order_file_name=4amber_${PDB}.order 4phenix_$PDB.updated.pdb


# echo
# '________________________________________________________Starting Model versus Data________________________________________________________'
# phenix.model_vs_data
# "${PDB}".updated_refine_001.pdb $PDB - sf.cif > "${PDB}"
# _model_v_data.txt
#
# echo
# '________________________________________________________Starting extract python script________________________________________________________'
# ~ / opt / anaconda3 / bin / python
# ~ / Fraser_Lab / phenix / parse_refine_log.py - log_file
# "${PDB}".updated_refine_001.log - PDB $PDB
# rm
# lig_RMSZ_updated.txt
# rm
# lig_RMSZ_pre_refine.txt
#
# echo
# '________________________________________________________Validating the Ligand from Original PDB________________________________________________________'
# if [[-e "${PDB}.ligands.cif"]]; then  # only running this if we have a ligand
#
# echo
# "$lig_name"
# # mmtbx.validate_ligands ${PDB}.pdb ${PDB}-sf.mtz ligand_code=$lig_name #"${PDB}_ligand" #prints out ADPs and occs + additional information
#
# echo
# '________________________________________________________Phenix PDB Interpretation_______________________________________________________'
# phenix.pdb_interpretation
# "${PDB}".updated.pdb
# "${PDB}".ligands.cif
# write_geo = True
#
# echo
# '________________________________________________________Elbow Refine_Geo_Display_______________________________________________________'
# elbow.refine_geo_display
# "${PDB}".updated.pdb.geo
# "$lig_name" >> lig_RMSZ_pre_refine.txt  # prints out deviations including ligand specific RMSD and RMSz values
# fi
#
# echo
# '________________________________________________________Validating the Ligand after Initial Refinement________________________________________________________'
# if [[-e "${PDB}.ligands.cif"]]; then  # only running this if we have a ligand
# # mmtbx.validate_ligands ${PDB}.updated_refine_001.pdb ${PDB}-sf.mtz ligand_code=$lig_name #"${PDB}_ligand" #prints out ADPs and occs + additional information
#
# echo
# '________________________________________________________Phenix PDB Interpretation_______________________________________________________'
# phenix.pdb_interpretation
# "${PDB}".updated_refine_001.pdb
# "${PDB}".ligands.cif
# write_geo = True
#
# echo
# '________________________________________________________Elbow Refine_Geo_Display_______________________________________________________'
# elbow.refine_geo_display
# "${PDB}".updated_refine_001.pdb.geo
# "$lig_name" >> lig_RMSZ_updated.txt
# # prints out deviations including ligand specific RMSD and RMSz values
# # calculate the energy difference of the ligand in the model and it relaxed RM1/AM1, but can be linked to 3rd party packages
# fi
#
# echo
# '________________________________________________________Validating Ligand Output_______________________________________________________'
# ~ / opt / anaconda3 / bin / python
# ~ / Fraser_Lab / phenix / ligand_geo_parser.py - pre_refine = lig_RMSZ_pre_refine.txt - post_refine = lig_RMSZ_updated.txt - PDB =$PDB
#
# #  echo '________________________________________________________Begin Ensemble Refinement_______________________________________________________'
# #  qsub /wynton/home/fraserlab/swankowicz/190419_Phenix_ensemble/grid_search_ens_refine.sh $PDB $3


def main():
    """ Parse command line args and run """

    parser = argparse.ArgumentParser(
        description="Runs initial refinement on a PDB file and prepares for ensemble refinement pipeline"
    )

    parser.add_argument("pdb_id", type=str, help="Accession ID for PDB of interest")

    args = parser.parse_args()

    pdb_id = args.pdb_id

    base_folder = '/Users/dmuir/kinome_ensemble_refinement/initial_egfr_refinements'

    # download_input_files(pdb_id, base_folder)
    #
    # # Source phenix environment
    # sp.call("source /Applications/phenix-1.16-3549/phenix_env.sh", shell=True)
    #
    # lig_name = parse_ligand_from_pdb(pdb_id, base_folder)

    os.chdir(base_folder + "/" + pdb_id)

    # print(lig_name)
    # print("\nStarting Phenix elbow...\n")
    # sp.call(f'phenix.elbow {pdb_id}.pdb --residue {lig_name} ', shell=True)
    #        # env={"PATH" :"/Applications/phenix-1.16-3549/build/bin/"}, shell=True)
    #
    # print("\nStarting Phenix cif as mtz\n")
    # sp.call(f'phenix.cif_as_mtz {pdb_id}-sf.cif --extend_flags --merge', shell=True)
    #
    # print("\nStarting Phenix Ready Set...\n")
    # sp.call(f'phenix.ready_set pdb_file_name={pdb_id}.pdb'
    #         f' cif_file_name=elbow.{lig_name}.{pdb_id}_pdb.001.cif >> readyset_output.txt', shell=True)

    run_refinement(pdb_id)


if __name__ == '__main__':
    main()
