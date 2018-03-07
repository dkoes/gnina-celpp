#!/usr/bin/env python

__author__ = 'j5wagner@ucsd.edu'


from d3r.celppade.custom_dock import Dock
import glob,re
import pybel, math

class gnina(Dock):
    """Abstract class defining methods for a custom docking solution
    for CELPP
    """
    Dock.SCI_PREPPED_LIG_SUFFIX = '_prepared.mol2'
    Dock.SCI_PREPPED_PROT_SUFFIX = '_prepared.pdb'

    def __init__(self,gnina_args=''):
        self.args = gnina_args
        Dock.__init__(self)

    def _ligand_length(self, lname):
        '''Return max distance between any two atoms of the ligand'''
        mol = pybel.readfile(os.path.splitext(lname)[1].lstrip('.'), lname).next()
        maxdist = 0
        for a in mol.atoms:
            for b in mol.atoms:
                sum = 0.0
                for i in xrange(3):
                    sum += (a.coords[i]-b.coords[i])**2
                dist = math.sqrt(sum)
                if dist > maxdist:
                    maxdist = dist
        return maxdist
        
    def ligand_technical_prep(self, sci_prepped_lig, targ_info_dict = {}):
        """
        'Technical preparation' is the step immediate preceding
        docking. During this step, you may perform any file
        conversions or processing that are specific to your docking
        program. Implementation of this function is optional.
        :param sci_prepped_lig: Scientifically prepared ligand file
        :param targ_info_dict: A dictionary of information about this target and the candidates chosen for docking.
        :returns: A list of result files to be copied into the
        subsequent docking folder. The base implementation merely
        returns the input string in a list (ie. [sci_prepped_lig]) 
        """
        return super(gnina,
                     self).ligand_technical_prep(sci_prepped_lig,
                                                 targ_info_dict = targ_info_dict)

    def receptor_technical_prep(self, 
                                sci_prepped_receptor, 
                                pocket_center, 
                                targ_info_dict = {}):
        """
        'Technical preparation' is the step immediately preceding
        docking. During this step, you may perform any file
        conversions or processing that are specific to your docking
        program. Implementation of this function is optional.
        :param sci_prepped_receptor: Scientifically prepared receptor file
        :param pocket_center: list of floats [x,y,z] of predicted pocket center
        :param targ_info_dict: A dictionary of information about this target and the candidates chosen for docking.
        :returns: A list of result files to be copied into the
        subsequent docking folder. This implementation merely
        returns the input string in a list (ie [sci_prepped_receptor])
        """
        threshold = 12

        #this is a bit ridiculous, have to figure out which receptor I'm
        #dealing with from the file name
        m = re.search(r'^([^-]+?)-',sci_prepped_receptor)
        prefix = m.group(1)
        #split receptor into receptor' and cognate ligand
        base = os.path.splitext(sci_prepped_receptor)[0]
        recname = base+'_rec.pdb'
                
        # punt on apo
        if 'lig_name' not in targ_info_dict[prefix][0]:
            os.system('cp %s %s'%(sci_prepped_receptor,recname))
            return [recname]
        
        lnames = set(targ_info_dict[prefix][0]['lig_name'])
        atomcnt = 0
        
        if lnames:
            pocketname = base+'_lig.pdb'
            rec = open(recname,'w')
            pocket = open(pocketname,'w')
            for line in open(sci_prepped_receptor):
                if (line.startswith('ATOM') or line.startswith('HETATM')) and line[17:20].strip() in lnames:
                    x = abs(float(line[30:38]) - pocket_center[0])
                    y = abs(float(line[38:46]) - pocket_center[1])
                    z = abs(float(line[46:54]) - pocket_center[2])
                    if x < threshold and y < threshold and z < threshold:
                        pocket.write(line)
                        atomcnt += 1
                else:
                    rec.write(line)
            rec.close()
            pocket.close()
            
        if atomcnt > 0:
            return [recname,pocketname]
        else:
            return [recname]


    def dock(self, 
             tech_prepped_lig_list, 
             tech_prepped_receptor_list, 
             output_receptor_pdb, 
             output_lig_mol, 
             targ_info_dict={}):
        """
        This function is the only one which the contestant MUST
        implement.  The dock() step runs the actual docking
        algorithm. Its first two arguments are the return values from
        the technical preparation functions for the ligand and
        receptor. These arguments are lists of file names (strings),
        which can be assumed to be in the current directory. 
        If prepare_ligand() and ligand_technical_prep() are not
        implemented by the contestant, tech_prepped_lig_list will
        contain a single string which names a SMILES file in the
        current directory.
        If receptor_scientific_prep() and receptor_technical_prep() are not
        implemented by the contestant, tech_prepped_receptor_list will
        contain a single string which names a PDB file in the current
        directory.
        The outputs from this step must be two files - a pdb with the
        filename specified in the output_receptor_pdb argument, and a
        mol with the filename specified in the output_ligand_mol
        argument.
        :param tech_prepped_lig_list: The list of file names resturned by ligand_technical_prep. These have been copied into the current directory.
        :param tech_prepped_receptor_list: The list of file names resturned by receptor_technical_prep. These have been copied into the current directory.
        :param output_receptor_pdb: The final receptor (after docking) must be converted to pdb format and have exactly this file name.
        :param output_lig mol: The final ligand (after docking) must be converted to mol format and have exactly this file name.
        :param targ_info_dict: A dictionary of information about this target and the candidates chosen for docking.
        :returns: True if docking is successful, False otherwise. Unless overwritten, this implementation always returns False
        """

        receptor = tech_prepped_receptor_list[0]
        ligand = tech_prepped_lig_list[0]
        
        ligandlen = self._ligand_length(ligand)
        # extract ligands from all identified similar receptors
        pocket_center = targ_info_dict['pocket_center']        
        output_lig = 'docked.sdf.gz'
        if len(tech_prepped_receptor_list) > 1: #receptor had ligand, use that
            pocketlig = tech_prepped_receptor_list[1]
            pocketlen = self._ligand_length(pocketlig)
            #increase autobox add if ligand is much bigger than pocket ligand
            add = max((ligandlen-pocketlen)/2.0, 4)
            gnina_command = 'gnina %s -r %s -l %s --autobox_ligand %s --autobox_add %f -o %s --seed 0 1> gnina.stdout 2> gnina.stderr' % (self.args, receptor, ligand, pocketlig, add, output_lig)
        else: #use box around ligand
            sz = max(16,ligandlen+4)
            gnina_command = 'gnina %s -r %s -l %s --center_x %f --center_y %f --center_z %f --size_x %f --size_y %f --size_z %f  -o %s --seed 0 1> gnina.stdout 2> gnina.stderr' % (self.args, receptor, ligand, pocket_center[0], pocket_center[1], pocket_center[2], sz, sz, sz, output_lig)

        print "Running: " + gnina_command
        os.system(gnina_command)
        os.system('babel %s %s' % (receptor, output_receptor_pdb))
        os.system('babel -l 1 %s %s' % (output_lig, output_lig_mol))
        return True




if ("__main__") == (__name__):
    import os
    import logging
    import shutil
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-l", "--ligsciprepdir", metavar="PATH", help = "PATH where we can find the scientific ligand prep output")
    parser.add_argument("-p", "--protsciprepdir", metavar="PATH", help = "PATH where we can find the scientific protein prep output")
    parser.add_argument("-o", "--outdir", metavar = "PATH", help = "PATH where we will put the docking output")
    parser.add_argument("--args",default='',help="additional arguments to provide to gnina")
    
    # Leave option for custom logging config here
    logger = logging.getLogger()
    logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.INFO )
    opt = parser.parse_args()
    lig_sci_prep_dir = opt.ligsciprepdir
    prot_sci_prep_dir = opt.protsciprepdir
    dock_dir = opt.outdir
    #running under this dir
    abs_running_dir = os.getcwd()
    log_file_path = os.path.join(abs_running_dir, 'final.log')
    log_file_dest = os.path.join(os.path.abspath(dock_dir), 'final.log')
    docker = gnina(opt.args)
    docker.run_dock(prot_sci_prep_dir,
                    lig_sci_prep_dir,
                    dock_dir)
    #move the final log file to the result dir
    shutil.move(log_file_path, log_file_dest)
