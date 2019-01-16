# XML equivalent
def configureSimulation(sim):
    import CompuCellSetup
    from XMLUtils import ElementCC3D
    
    CompuCell3DElmnt=ElementCC3D("CompuCell3D",{"Revision":"20170130","Version":"3.7.7"})
    PottsElmnt=CompuCell3DElmnt.ElementCC3D("Potts")
    PottsElmnt.ElementCC3D("Dimensions",{"x":"4","y":"20","z":"1"})
    PottsElmnt.ElementCC3D("Steps",{},"28801")
    PottsElmnt.ElementCC3D("Temperature",{},"10.0")
    PottsElmnt.ElementCC3D("NeighborOrder",{},"1")
    MetadataElmnt=CompuCell3DElmnt.ElementCC3D("Metadata")
    MetadataElmnt.ElementCC3D("DebugOutputFrequency",{},"500")
    MetadataElmnt.ElementCC3D("NumberOfProcessors",{},"2")
    PluginElmnt=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"CellType"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"0","TypeName":"Medium"})
    PluginElmnt.ElementCC3D("CellType",{"Freeze":"","TypeId":"1","TypeName":"Blood"})
    PluginElmnt.ElementCC3D("CellType",{"Freeze":"","TypeId":"2","TypeName":"Hepatocyte"})
    PluginElmnt_1=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"Volume"})
    PluginElmnt_1.ElementCC3D("VolumeEnergyParameters",{"CellType":"Blood","LambdaVolume":"2.0","TargetVolume":"200"})
    PluginElmnt_1.ElementCC3D("VolumeEnergyParameters",{"CellType":"Hepatocyte","LambdaVolume":"2.0","TargetVolume":"400"})

    CompuCellSetup.setSimulationXMLDescription(CompuCell3DElmnt)    

# Main Python script    
import sys
from os import environ
from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])

import CompuCellSetup

sim,simthread = CompuCellSetup.getCoreSimulationObjects()
            
# add extra attributes here

# Call XML replacement function            
configureSimulation(sim)
            
CompuCellSetup.initializeSimulationObjects(sim,simthread)          

# Definitions of additional Python-managed fields go here
        
#Add Python steppables here
steppableRegistry=CompuCellSetup.getSteppableRegistry()
        
from multi_simple_liverSteppables import multi_simple_liverSteppable
steppableInstance=multi_simple_liverSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(steppableInstance)
   
CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
