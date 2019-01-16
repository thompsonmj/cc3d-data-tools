#! /usr/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anshaikh@iu.edu
#SBATCH -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/stdout.txt
#SBATCH --nodes=8 
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_10/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_10/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_10/output -f 50
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_20/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_20/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_20/output -f 50
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_30/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_30/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_30/output -f 50
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_40/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_40/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_40/output -f 50
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_50/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_50/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulationExtInp/OutputEI/Sim_50/output -f 50
