#!/bin/bash
# Clear previous results so we don't get mixed up
rm -r 1-get_challenge_data
rm -r 2-protein_prep
rm -r 3-ligand_prep
rm -r 4-docking
rm -r 4-docking_cnn
rm -r 5-pack_docking_results
rm -r 5-pack_docking_results_cnn

source /home/dkoes/.bashrc

# Run CELPP for this week
mkdir 1-get_challenge_data
getchallengedata.py --unpackdir 1-get_challenge_data -f ../ftp_config

mkdir 2-protein_prep
python ../gnina-celpp/gnina-celpp_protein_prep.py --challengedata 1-get_challenge_data/ --prepdir 2-protein_prep

mkdir 3-ligand_prep
python ../gnina-celpp/gnina-celpp_ligand_prep.py --challengedata 1-get_challenge_data/ --prepdir 3-ligand_prep

mkdir 4-docking
python ../gnina-celpp/gnina-celpp_dock.py --protsciprepdir 2-protein_prep --ligsciprepdir 3-ligand_prep --outdir  4-docking

mkdir 5-pack_docking_results
packdockingresults.py --dockdir 4-docking --packdir 5-pack_docking_results -f ../ftp_config --challengedata 1-get_challenge_data/


mkdir 4-docking_cnn
python ../gnina-celpp/gnina-celpp_dock.py --args "--cnn_scoring  --minimize_early_term --gpu" --protsciprepdir 2-protein_prep --ligsciprepdir 3-ligand_prep --outdir  4-docking_cnn

mkdir 5-pack_docking_results_cnn
packdockingresults.py --dockdir 4-docking_cnn --packdir 5-pack_docking_results_cnn -f ../ftp_config2 --challengedata 1-get_challenge_data/
