#! /bin/sh -l
#PBS -N CompuCell3D
#PBS -q dumulis 
#PBS -M anshaikh@iu.edu
#PBS -m abe
#PBS -l nodes=2:ppn=2
 
qsub /home/thompsj/CompuCell3D/CC3D_v377/CLI/CompuCell3DCLI.py -i /home/thompsmj/CompuCell3D/TestSimulationExtInp/OutputBrown/Sim_10/TestSimulation1.cc3d -o /home/thompsmj/CompuCell3D/TestSimulationExtInp/OutputBrown/Sim_10/output -f 50

qsub /home/thompsj/CompuCell3D/CC3D_v377/CLI/CompuCell3DCLI.py -i /home/thompsmj/CompuCell3D/TestSimulationExtInp/OutputBrown/Sim_20/TestSimulation1.cc3d -o /home/thompsmj/CompuCell3D/TestSimulationExtInp/OutputBrown/Sim_20/output -f 50

