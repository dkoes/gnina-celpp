#!/usr/bin/env python


__author__ = 'dkoes@pitt.edu'

from d3r.celppade.custom_protein_prep import ProteinPrep
import os, subprocess


class dkoes_protprep(ProteinPrep):
    """Abstract class defining methods for a custom docking solution
    for CELPP
    """
    ProteinPrep.OUTPUT_PROTEIN_SUFFIX = '.pdb'
    def receptor_scientific_prep(self, 
                                 protein_file, 
                                 prepared_protein_file, 
                                 targ_info_dict={}):
        """
        Protein 'scientific preparation' is the process of generating
        a dockable representation of the candidate protein from a
        single-chain PDB file.
        :param protein_file: PDB file containing candidate protein.  
        :param prepared_protein_file: The result of preparation should have this file name.  
        :param targ_info_dict: A dictionary of information about this target and the candidates chosen for docking.  
        :returns: True if preparation was successful. False otherwise.
        """
        #remove alt locations and water
        subprocess.check_call('pdb4amber -d %s -o %s'%(protein_file,prepared_protein_file),shell=True)
        return True

 

    
if ("__main__") == (__name__):
    import logging
    import os
    import shutil
    from argparse import ArgumentParser
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

    prot_prepper = dkoes_protprep()
    prot_prepper.run_scientific_protein_prep(challenge_data_path, pdb_location, prep_result_path)

    #move the final log file to the result dir
    shutil.move(log_file_path, log_file_dest)
