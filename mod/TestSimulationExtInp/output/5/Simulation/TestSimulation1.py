
def configureSimulation(sim):
    import CompuCellSetup
    from XMLUtils import ElementCC3D
    from Parameters import *
    CompuCell3DElmnt=ElementCC3D("CompuCell3D",{"Revision":"20170723","Version":"3.7.6"})
    
    PottsElmnt=CompuCell3DElmnt.ElementCC3D("Potts")
    
    # Basic properties of CPM (GGH) algorithm
    PottsElmnt.ElementCC3D("Dimensions",{"x":"100","y":"100","z":"1"})
    PottsElmnt.ElementCC3D("Steps",{}, totalNumberOfMCS)
    PottsElmnt.ElementCC3D("Temperature",{},temperature)
    PottsElmnt.ElementCC3D("NeighborOrder",{},neighborOrder)
    
    PluginElmnt=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"CellType"})
    
    # Listing all cell types in the simulation
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"0","TypeName":"Medium"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"1","TypeName":"CellA"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"2","TypeName":"CellB"})
    
    PluginElmnt_1=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"CenterOfMass"})
    
    # Module tracking center of mass of each cell
    
    PluginElmnt_2=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"NeighborTracker"})
    
    # Module tracking neighboring cells of each cell
    
    PluginElmnt_3=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"Contact"})
    # Specification of adhesion energies
    PluginElmnt_3.ElementCC3D("Energy",{"Type1":"Medium","Type2":"Medium"}, medium_medium)
    PluginElmnt_3.ElementCC3D("Energy",{"Type1":"Medium","Type2":"CellA"}, medium_cellA)
    PluginElmnt_3.ElementCC3D("Energy",{"Type1":"Medium","Type2":"CellB"}, medium_cellB)
    PluginElmnt_3.ElementCC3D("Energy",{"Type1":"CellA","Type2":"CellA"}, cellA_cellA)
    PluginElmnt_3.ElementCC3D("Energy",{"Type1":"CellA","Type2":"CellB"}, cellA_cellB)
    PluginElmnt_3.ElementCC3D("Energy",{"Type1":"CellB","Type2":"CellB"}, cellB_cellB)
    PluginElmnt_3.ElementCC3D("NeighborOrder",{},"1")
    
    SteppableElmnt=CompuCell3DElmnt.ElementCC3D("Steppable",{"Type":"UniformInitializer"})
    
    # Initial layout of cells in the form of rectangular slab
    RegionElmnt=SteppableElmnt.ElementCC3D("Region")
    RegionElmnt.ElementCC3D("BoxMin",{"x":"20","y":"20","z":"0"})
    RegionElmnt.ElementCC3D("BoxMax",{"x":"80","y":"80","z":"1"})
    RegionElmnt.ElementCC3D("Gap",{},"0")
    RegionElmnt.ElementCC3D("Width",{},"5")
    RegionElmnt.ElementCC3D("Types",{},"CellA,CellB")

    CompuCellSetup.setSimulationXMLDescription(CompuCell3DElmnt)    
    


            
    CompuCellSetup.setSimulationXMLDescription(CompuCell3DElmnt)
            
import sys
from os import environ
from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])


import CompuCellSetup


sim,simthread = CompuCellSetup.getCoreSimulationObjects()
        
configureSimulation(sim)            
            
# add extra attributes here
        
CompuCellSetup.initializeSimulationObjects(sim,simthread)
# Definitions of additional Python-managed fields go here
        
#Add Python steppables here
steppableRegistry=CompuCellSetup.getSteppableRegistry()
        
from TestSimulation1Steppables import TestSimulation1Steppable
steppableInstance=TestSimulation1Steppable(sim,_frequency=1)
steppableRegistry.registerSteppable(steppableInstance)
        
CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        