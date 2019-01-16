#! /usr/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anshaikh@iu.edu
#SBATCH -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/stdout.txt
#SBATCH --nodes=8 
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/0/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/0/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/0/output -f 20
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/1/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/1/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/1/output -f 20
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/2/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/2/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/2/output -f 20
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/3/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/3/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/3/output -f 20
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/4/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/4/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/4/output -f 20
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/5/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/5/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/5/output -f 20
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/6/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/6/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/6/output -f 20
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/7/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/7/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/7/output -f 20
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/8/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/8/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/8/output -f 20
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/9/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/9/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/9/output -f 20
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/10/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/10/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/10/output -f 20
srun -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/11/output.txt CC3D_CLI -i /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/11/TestSimulation1.cc3d -o /Users/anwar/CC3D/Workshop/CC3D/TestSimulation1/output/11/output -f 20
