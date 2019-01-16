sim=/depot/dumulis/cc3d/apps/cc3d-data-integration-and-optimization/mod/multi_simple_liver/multi_simple_liver.cc3d
out=/depot/dumulis/cc3d/apps/cc3d-data-integration-and-optimization/mod/multi_simple_liver/testout
par=/depot/dumulis/cc3d/apps/cc3d-data-integration-and-optimization/mod/multi_simple_liver/par_sgL2.csv

python ../../dep/CompuCell3D/CompuCell3D/CLI/CompuCell3DCLI.py -i $sim -o $out -p $par -f 100
#python ../../../cc3d_v377/CLI/CompuCell3DCLI.py -i $sim -o $out -p $par -f 100
#python ../../../old/dep/cc3d_v377/CLI/CompuCell3DCLI.py -i $sim -o $out -p $par -f 100
#python /home/thompsmj/CC3D/CC3D_v377/CLI/CompuCell3DCLI.py -i $sim -o $out -p $par -f 100
