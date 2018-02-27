#!/usr/bin/env python

__author__ = 'dkoes@pitt.edu'

from d3r.celppade.custom_ligand_prep import LigandPrep
import os
import glob, re

rdkit_smiles_to_3d_sdf_text = '''
import rdkit.Chem
import rdkit.Chem.AllChem
import sys
if not(len(sys.argv)) == 3:
    print "python smiles2Mol.py inputSmiles outputSdf"
    sys.exit()
smiles = open(sys.argv[1]).read().strip()
mol = rdkit.Chem.MolFromSmiles(smiles)
molH = rdkit.Chem.AddHs(mol)
rdkit.Chem.AllChem.SanitizeMol(molH)
rdkit.Chem.AllChem.EmbedMolecule(molH)
rdkit.Chem.AllChem.UFFOptimizeMolecule(molH)
w = rdkit.Chem.SDWriter(sys.argv[2])
w.write(molH)
w.close()
'''




class dkoes_ligprep(LigandPrep):
            
    """Class defining methods for a custom ligand docking solution
    for CELPP
    """
    LigandPrep.OUTPUT_LIG_SUFFIX = '.mol2'
    def ligand_scientific_prep(self, 
                               lig_smi_file, 
                               out_lig_file, 
                               targ_info_dict={}):
        """
        Ligand 'scientific preparation' is the process of generating a
        dockable representation of the target ligand from its SMILES
        string.
        :param lig_smi_file: File containing SMILES for target ligand.  
        :param out_lig_file: The result of preparation should have this file name.  
        :param targ_info_dict: A dictionary of information about this target and the candidates chosen for docking.  
        :returns: True if preparation was successful. False otherwise.
        """

        with open('rdkit_smiles_to_3d_sdf.py','wb') as of:
            of.write(rdkit_smiles_to_3d_sdf_text)

        os.system('python rdkit_smiles_to_3d_sdf.py ' + lig_smi_file + 
                  ' ligand.sdf 1> rdkit_smiles_to_3d_sdf.stdout 2>' +
                  ' rdkit_smiles_to_3d_sdf.stderr')
        os.system('babel -isdf ligand.sdf -omol2 ligand.mol2'
                  + ' 1> lig_sdf_to_mol2.stdout 2> lig_sdf_to_mol2.stderr')
        os.system('cp ligand.mol2 ' + out_lig_file)

        return True








if ("__main__") == (__name__):
    from argparse import ArgumentParser
    import os
    import logging 
    import shutil
    parser = ArgumentParser()
    parser.add_argument("-p", "--pdbdb", metavar = "PATH", help = "PDB DATABANK which we will dock into")
    parser.add_argument("-c", "--challengedata", metavar="PATH", help = "PATH to the unpacked challenge data package")
    parser.add_argument("-o", "--prepdir", metavar = "PATH", help = "PATH to the output directory")
    logger = logging.getLogger()
    logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level = logging.INFO )
    opt = parser.parse_args()
    pdb_location = opt.pdbdb
    challenge_data_path = opt.challengedata
    prep_result_path = opt.prepdir

    #running under this dir
    abs_running_dir = os.getcwd()
    log_file_path = os.path.join(abs_running_dir, 'final.log')
    log_file_dest = os.path.join(os.path.abspath(prep_result_path), 'final.log')

    lig_prepper =  dkoes_ligprep()
    lig_prepper.run_scientific_ligand_prep(challenge_data_path, pdb_location, prep_result_path)

    #move the final log file to the result dir
    shutil.move(log_file_path, log_file_dest)

