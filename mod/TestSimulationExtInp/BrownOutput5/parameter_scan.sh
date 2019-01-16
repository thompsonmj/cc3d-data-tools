#! /bin/sh -l
#PBS -q dumulis
#PBS -M thompsmj@purdue.edu
#PBS -l nodes=4:ppn=4
#PBS -m abe
#PBS -N CompuCell3Dqsub /home/thompsmj/CompuCell3D/CC3D_v377/CLI/CompuCell3DCLI.py -i /home/thompsmj/CompuCell3D/TestSimulationExtInp/BrownOutput5/Sim_10/TestSimulation1.cc3d -o /home/thompsmj/CompuCell3D/TestSimulationExtInp/BrownOutput5/Sim_10/output -f 50
qsub /home/thompsmj/CompuCell3D/CC3D_v377/CLI/CompuCell3DCLI.py -i /home/thompsmj/CompuCell3D/TestSimulationExtInp/BrownOutput5/Sim_20/TestSimulation1.cc3d -o /home/thompsmj/CompuCell3D/TestSimulationExtInp/BrownOutput5/Sim_20/output -f 50
qsub /home/thompsmj/CompuCell3D/CC3D_v377/CLI/CompuCell3DCLI.py -i /home/thompsmj/CompuCell3D/TestSimulationExtInp/BrownOutput5/Sim_30/TestSimulation1.cc3d -o /home/thompsmj/CompuCell3D/TestSimulationExtInp/BrownOutput5/Sim_30/output -f 50
qsub /home/thompsmj/CompuCell3D/CC3D_v377/CLI/CompuCell3DCLI.py -i /home/thompsmj/CompuCell3D/TestSimulationExtInp/BrownOutput5/Sim_40/TestSimulation1.cc3d -o /home/thompsmj/CompuCell3D/TestSimulationExtInp/BrownOutput5/Sim_40/output -f 50
qsub /home/thompsmj/CompuCell3D/CC3D_v377/CLI/CompuCell3DCLI.py -i /home/thompsmj/CompuCell3D/TestSimulationExtInp/BrownOutput5/Sim_50/TestSimulation1.cc3d -o /home/thompsmj/CompuCell3D/TestSimulationExtInp/BrownOutput5/Sim_50/output -f 50
