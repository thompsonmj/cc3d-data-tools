from PlayerPython import * 
import CompuCellSetup
from PySteppables import *
import CompuCell
import sys
import math
import os

###import multi_simple_liver_parameters as p
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# total simulation time 
totalSimTimeHr = 1 # 8 hours default, in vitro data goes out to 24 hours

# number of iterations per mcs on the CC3D's model ODEs
cc3d_iterODE = 1  # best = 5, Default = 1. If 5, then CC3D will divide the transfer ODE calc's into 5 substeps

# conversion factor for CC3D mcs to simulated time
# note that if this value is changed the total number of mcs needs to be changed as well
cc3d_mcs2second = 1.0 # Best = 0.5, Default = 1.0 simulated second per mcs (second/mcs), if 0.5 then 2 mcs per simulated second


# File: multi_simple_liver_paramters.py
# Date:  1 Oct 2017
# Created By: J Sluka
# Comments:
# 
# C:\Users\James Sluka\Desktop\Small Code Projects\multi_simple_liver\
# 
#######################################################################################################################
# Name = value	# Description	Unit	Ref.	Ontology Term URI
# Subcellular Metabolism of APAP -> APAPG + APAPS #####################################################################
APAP=0.0	# Acetaminophen  gram    					
NAPQI=0.0	# N-acetyl para amino quinone   gram					
GSH=10.0	# Glutathione in hepatocytes    mM					
NAPQIGSH=0.0	# APAP+NAPQI conjugate  gram					
APAPconj_Glu=0.0	# APAP-Glucuronide conjugate    gram					
APAPconj_Sul=0.0	# APAP-Sulfate conjugate    gram					
X1=0.0	# 	
##### N E W ##################### Jan 2018
#  hepatocyte volume =8x10^-15L  (20um)^3
#  sinusoid segment volume is 1.3x10^-15L (14% of hep + sinuoid)
#  ratio = 8/1.3 = 6.15, which is the volume ratio of heps to sinusoid
hepatoVolume = 8.0e-15 # Liter
bloodVolume  = 1.3e-15 # Liter, this is the blood volume in a hep long sinusoid,
                       # based on ~15% of liver volume being sinusoid
##### N E W ##################### Jan 2018
# [APAP] 0.0	# 					
# [NAPQI] 0.0	# 					
# [GSH] 10.0	# 					
# [NAPQIGSH] 0.0	# 					
# [APAPconj_Glu] 0.0	# 					
# [APAPconj_Sul] 0.0	# 					
# [X1] 0.0	# 					
# compartment 1.0	# 					
# Vmax_2E1_APAP 2e-05	# 					
# Km_2E1_APAP 1.29	# 					
# kNapqiGsh 0.1	# 					
# kGsh 0.0001	# 					
# GSHmax 10.0	# 					
# Vmax_PhaseIIEnzGlu_APAP 0.001	# 					
# Km_PhaseIIEnzGlu_APAP 1.0	# 					
# Vmax_PhaseIIEnzSul_APAP 0.000175	# 					
# Km_PhaseIIEnzSul_APAP 0.2	# 					
# init(Vmax_2E1_APAP) 2e-05	# 					
# init(Km_2E1_APAP) 1.29	# 					
# init(kNapqiGsh) 0.1	# 					
# init(kGsh) 0.0001	# 					
# init(GSHmax) 10.0	# 					
# init(Vmax_PhaseIIEnzGlu_APAP) 0.001	# 					
# init(Km_PhaseIIEnzGlu_APAP) 1.0	# 					
# init(Vmax_PhaseIIEnzSul_APAP) 0.000175	# 					
# init(Km_PhaseIIEnzSul_APAP) 0.2	# 					
# J0 0.0	# 					
# J1 0.0	# 					
# J2 0.0	# 					
# J3 0.0	# 					
# J4 0.0	# 					
# init([APAP]) 0.0	# 					
# init([NAPQI]) 0.0	# 					
# init([GSH]) 10.0	# 					
# init([NAPQIGSH]) 0.0	# 					
# init([APAPconj_Glu]) 0.0	# 					
# init([APAPconj_Sul]) 0.0	# 					
# init(APAP) 0.0	# 					
# init(NAPQI) 0.0	# 					
# init(GSH) 10.0	# 					
# init(NAPQIGSH) 0.0	# 					
# init(APAPconj_Glu) 0.0	# 					
# init(APAPconj_Sul) 0.0	# 					
# APAP' 0.0	# 					
# NAPQI' 0.0	# 					
# GSH' 0.0	# 					
# NAPQIGSH' 0.0	# 					
# APAPconj_Glu' 0.0	# 					
# APAPconj_Sul' 0.0	# 					
# 
####################################################################################################################### 
# CC3D Model ##########################################################################################################
# passive rate constant blood-hep
# hep-hep
# active vmax
# active km
####################################################################################################################### 
# ALL THE PARAMETERS IN THE LIVER I MODEL with their value in the HMP sim6 set (sum of RMSE 1.791) ###################
# Parameter	  Value	# Unit
pbpk_Fup	= 0.800	 # unitless
pbpk_FupG	= 1.000	 # unitless
pbpk_FupS	= 1.000	 # unitless
pbpk_Kl2p	= 1.000	 # unitless  Kliver2plasma
pbpk_Kl2pG	= 1.000	 # unitless
pbpk_Kl2pS	= 1.000	 # unitless
pbpk_Kk2p	= 1.000	 # unitless  Kkidney2plasma
pbpk_Kk2pG	= 1.000	 # unitless
pbpk_Kk2pS	= 1.000	 # unitless
pbpk_Kr2p	= 1.700	 # unitless  KRest2plasma
pbpk_Kr2pG	= 0.250	 # unitless
pbpk_Kr2pS	= 0.100	 # unitless
pbpk_Qgfr	= 0.714	 # L/hr
pbpk_QgfrG	= 7.86	 # L/hr   7.86  12 15 18 22##############
pbpk_QgfrS	= 9.96	 # L/hr   9.96  15 20 22  ##############
pbpk_Rb2p	= 1.090	 # unitless
pbpk_Rb2pG	= 0.550	 # unitless
pbpk_Rb2pS	= 0.550	 # unitless
pbpk_bw         = 70.0   # Kg
pbpk_dose	= 1.400	 # g
pbpk_hemat	= 0.450	 # unitless
pbpk_kGutabs	= 1.400	 # 1/hr		
pbpk_CLmetabolism=0.0    # 1/hr0.000600

cc3d_Km_AT_APAP	= 0.010	 # mmol/L  0.010
cc3d_Vmax_AT_APAP=0.010	 # mmol/L/s  0.010  0.005 0.00010 0.00030 #######1/26/18#######
cc3d_k_AT_APAPG	= 0.000450 # 1/s
cc3d_k_AT_APAPS	= 0.001900 # 1/s
cc3d_k_PD_H2H	= 0.00100  # 1/s
cc3d_k_PD_H2R	= 0.00100  # 1/s
cc3d_k_PD_H2S	= 0.00100  # 1/s
cc3d_k_PD_R2H	= 0.00100  # 1/s
cc3d_k_PD_R2R	= 0.00100  # 1/s
cc3d_k_PD_R2S	= 0.00100  # 1/s
cc3d_k_PD_S2H	= 0.00100  # 1/s
cc3d_k_PD_S2R	= 0.00120  # 1/s
cc3d_k_PD_S2S	= 0.010	   # 1/s

sc_Km_2E1_APAP	= 1.290	   # mmol/L
sc_Km_GLUC	= 1.500	   # mmol/L
sc_Km_SULF	= 0.300	   # mmol/L
sc_Vmax_2E1_APAP= 0.000020 # mmol/L/s   0.000020
sc_Vmax_GLUC	= 0.000900 # mmol/L/s
sc_Vmax_SULF	= 0.000150 # mmol/L/s
sc_kGsh	        = 0.000100 # mmol/L/s
sc_kNapqiGsh	= 0.100	   # mmol/L/s
####### RESULTS FOR THIS PARAMETER SET
# RMSEsum  =  1.791090615
# RMSEapap  =  1.032493107
# RMSEapapg  =  0.351259765
# RMSEapaps  =  0.407337742
# CmaxA  =  15.21936172
# CmaxG  =  10.85402564
# CmaxS  =  4.488459375
# TmaxA  =  1.308333333
# TmaxG  =  3.475
# TmaxS  =  2.425
# AUCA  =  67.5451445
# AUCG  =  65.31191228
# AUCS  =  27.69530186
# metabRatio  =  90.15959711
####################################################################################################################### 
# APAP PBPK Model #####################################################################################################
# CArt 0.0	# 					
# CGut 0.0	# 					
# AGutlumen 0.00925925925926	# 					
# CLung 0.0	# 					
# CVen 0.0	# 					
# CRest 0.0	# 					
# CLiver 0.0	# 					
# CMetabolized 0.0	# 					
# CKidney 0.0	# 					
# CTubules 0.0	# 					
# [CArt] 0.0	# 					
# [CGut] 0.0	# 					
# [AGutlumen] 0.00925925925926	# NOTE THAT [AGutlumen] is reported to be the same as AGutlumen, and both are concentrations		
# [CLung] 0.0	# 					
# [CVen] 0.0	# 					
# [CRest] 0.0	# 					
# [CLiver] 0.0	# 					
# [CMetabolized] 0.0	# 					
# [CKidney] 0.0	# 					
# [CTubules] 0.0	# 					
# VArt 1.4994	# 					
# VGut 1.1046	# 					
# VGutLumen 1.0	# 					
# VLung 0.5082	# 					
# VVen 3.4104	# 					
# VRest 33.4698	# 					
# VLiver=1.7136*0.14	# overriding the calculation of the liver, set to 14% of the liver volume to represent just the blood volume				 					
# VKidney 0.294	# 					
# VKidneyTubules 1.0	# 					
# APAP_Dose_grams 1.4	# 					
# APAP_MW 151.2	# 					
# APAP_Dose 0.00925925925926	# 					
# BW 70.0	# 					
# BW_ref 1.0	# 					
# Cardiac_flow_ref 15.0	# 					
# QCardiac 363.006823874	# 					
# QGut_fraction_QCardiac 0.205	# 					
# QGut 74.4163988942	# 					
# QLiver_fraction_QCardiac 0.0535	# 					
# QLiver 19.4208650773	# 					
# QKidney_fraction_QCardiac 0.2214	# 					
# QKidney 80.3697108057	# 					
# QRest 188.799849097	# 					
# QGFR_ref 0.039	# 					
# Qgfr 0.943817742072	# 					
# kGutabs 1.5	# 					
# CLmetabolism 0.0	# 					
# Fraction_unbound_plasma 0.8	# 					
# Ratioblood2plasma 1.09	# 					
# Kliver2plasma =1.0	# 					
# Kkidney2plasma=1.0	# 					
# KRest2plasma  =1.6	# 					
# VTotal_ref 0.6	# 					
# VTotal 42.0	# 					
# VLiver_fraction_VTotal 0.0408
# VArt_fraction_VTotal 0.0357	# 					
# VLung_fraction_VTotal 0.0121	# 					
# VVen_fraction_VTotal 0.0812	# 					
# VGut_fraction_VTotal 0.0263	# 					
# VKidney_fraction_VTotal 0.007	# 					
# init(APAP_Dose_grams) 1.4	# 					
# init(APAP_MW) 151.2	# 					
# init(APAP_Dose) 0.00925925925926	# 					
# init(BW) 70.0	# 					
# init(BW_ref) 1.0	# 					
# init(Cardiac_flow_ref) 15.0	# 					
# init(QCardiac) 363.006823874	# 					
# init(QGut_fraction_QCardiac) 0.205	# 					
# init(QGut) 74.4163988942	# 					
# init(QLiver_fraction_QCardiac) 0.0535	# 					
# init(QLiver) 19.4208650773	# 					
# init(QKidney_fraction_QCardiac) 0.2214	# 					
# init(QKidney) 80.3697108057	# 					
# init(QRest) 188.799849097	# 					
# init(QGFR_ref) 0.039	# 					
# init(Qgfr) 0.943817742072	# 					
# init(kGutabs) 1.5	# 					
# init(CLmetabolism) 9.5	# 					
# init(Fraction_unbound_plasma) 0.8	# 					
# init(Ratioblood2plasma) 1.09	# 					
# init(Kliver2plasma) 1.0	# 					
# init(Kkidney2plasma) 1.0	# 					
# init(KRest2plasma) 1.6	# 					
# init(VTotal_ref) 0.6	# 					
# init(VTotal) 42.0	# 					
# init(VLiver_fraction_VTotal) 0.0408	# 					
# init(VArt_fraction_VTotal) 0.0357	# 					
# init(VLung_fraction_VTotal) 0.0121	# 					
# init(VVen_fraction_VTotal) 0.0812	# 					
# init(VGut_fraction_VTotal) 0.0263	# 					
# init(VKidney_fraction_VTotal) 0.007	# 					
# J1 0.0	# 					
# J2 0.0138888888889	# 					
# J3 0.0	# 					
# J4 0.0	# 					
# J5 0.0	# 					
# J6 0.0	# 					
# J7 0.0	# 					
# J8 0.0	# 					
# J9 0.0	# 					
# J10 0.0	# 					
# J11 0.0	# 					
# J12 0.0	# 					
# J13 0.0	# 					
# init([CArt]) 0.0	# 					
# init([CGut]) 0.0	# 					
# init([AGutlumen]) 0.00925925925926	# 					
# init([CLung]) 0.0	# 					
# init([CVen]) 0.0	# 					
# init([CRest]) 0.0	# 					
# init([CLiver]) 0.0	# 					
# init([CMetabolized]) 0.0	# 					
# init([CKidney]) 0.0	# 					
# init([CTubules]) 0.0	# 					
# init(CArt) 0.0	# 					
# init(CGut) 0.0	# 					
# init(AGutlumen) 0.00925925925926	# 					
# init(CLung) 0.0	# 					
# init(CVen) 0.0	# 					
# init(CRest) 0.0	# 					
# init(CLiver) 0.0	# 					
# init(CMetabolized) 0.0	# 					
# init(CKidney) 0.0	# 					
# init(CTubules) 0.0	# 					
# CArt' 0.0	# 					
# CGut' 0.0138888888889	# 					
# AGutlumen' -0.0138888888889	# 					
# CLung' 0.0	# 					
# CVen' 0.0	# 					
# CRest' 0.0	# 					
# CLiver' 0.0	# 					
# CMetabolized' 0.0	# 					
# CKidney' 0.0	# 					
# CTubules' 0.0	# 					
 
 
# # # ## ## ##
# # # pbpk_Fup	= 0.7	 # unitless
# # # pbpk_Kk2p	= 1.	 # unitless  Kkidney2plasma 
# # # pbpk_Kr2p	= 1.8	 # unitless  KRest2plasma

# # # pbpk_Qgfr	= 0.61	 # L/hr
# # # pbpk_QgfrG	= 7.86	 # L/hr   7.86  12 15 18 22##############
# # # pbpk_QgfrS	= 9.96	 # L/hr   9.96  15 20 22  ##############

# # # pbpk_kGutabs	= 1.7	 # 1/hr	

cc3d_k_AT_APAPG = 0.00005 # 1/s  0.000450  1/30
cc3d_k_AT_APAPS = 0.00040 # 1/s  0.001900

cc3d_Km_AT_APAP   = 0.00200	 # mmol/L  0.010
cc3d_Vmax_AT_APAP = 0.00060	 # mmol/L/s  0.010  0.005 0.00010 0.00030 #######1/26/18#######

# cc3d_k_PD_H2S	= 0.000050       # 1/s   0.00100
# cc3d_k_PD_S2H	= cc3d_k_PD_H2S  # 1/s  0.00100

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class multi_simple_liverSteppable(SteppableBasePy):    
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.scalarCL_APAP=self.createScalarFieldCellLevelPy("APAP")            
        self.scalarCL_APAPG=self.createScalarFieldCellLevelPy("APAPG")            
        self.scalarCL_APAPS=self.createScalarFieldCellLevelPy("APAPS")            
            
    def start(self):
        # output file
        self.fileHandleLog, fullFileName = self.openFileInSimulationOutputDirectory("data_log.txt", "w")
        titleLine = "mcs\ttime(hr)"
        titleLine += "\thep1.dict['APAP']\thep1.dict['APAPG']\thep1.dict['APAPS']"
        titleLine += "\tblood1.dict['APAP']\tblood1.dict['APAPG']\tblood1.dict['APAPS']"
        titleLine += "\tPBPKspeciesDictAPAP['CLiver']\tPBPKspeciesDictAPAPG['CLiver']"
        titleLine += "\tPBPKspeciesDictAPAPS['CLiver']\tPBPKspeciesDictAPAP['CVen']"
        titleLine += "\tPBPKspeciesDictAPAPG['CVen']\tPBPKspeciesDictAPAPS['CVen']"
        titleLine += "\tPBPKspeciesDictAPAP['CTubules']\tPBPKspeciesDictAPAPG['CTubules']"
        titleLine += "\tPBPKspeciesDictAPAPS['CTubules']\tSUBCELLspeciesDict['APAP']"
        titleLine += "\tSUBCELLspeciesDict['APAPconj_Glu']\tSUBCELLspeciesDict['APAPconj_Sul']"
        titleLine += "\tSUBCELLspeciesDict['NAPQI']\tSUBCELLspeciesDict['NAPQIGSH']"
        titleLine += "\tSUBCELLspeciesDict['GSH']\n"
        self.fileHandleLog.write(titleLine)
            
        # make a hepatocyte and a blood portion "cell"
        # relative areas (volumes) of 15:85, though it  doesn't matter in this model
        self.cellField[ 0:4, 0:15,0:1] = self.newCell(self.HEPATOCYTE) 
        self.cellField[ 0:4,15:19,0:1] = self.newCell(self.BLOOD) 
      
        for cell in self.cellList:
            # access/modification of a dictionary attached to cell
            cell.dict['APAP']  = 0.0
            cell.dict['APAPG'] = 0.0
            cell.dict['APAPS'] = 0.0
            self.scalarCL_APAP[cell] =cell.dict['APAP']
            self.scalarCL_APAPG[cell]=cell.dict['APAPG'] 
            self.scalarCL_APAPS[cell]=cell.dict['APAPS']
                    
        #################################################
        #  time step, second/mcs
        self.mcs2second = cc3d_mcs2second # 1.0 simulated second per mcs (second/mcs)

        # "cc3d_iterODE" (in the paramter file) defins how man  steps while doing the CC3D diffusion ODEs, 
        # "2" means the ODEs are calculated in two steps, each of 0.5*mcs2seconds
        
        # calculate the number of MCS based on mcs2seconds value and total simulation time totalSimTimeHr.
        # (Note that "cc3d_iterODE" doesn't affect this
        self.setMaxMCS(int(3600*totalSimTimeHr/self.mcs2second + 1))   # new, Jan 2018

        # global parameters for the SBML solver
       #SBMLoptions = {'relative': 1e-10, 'absolute': 1e-12, 'steps': 10, 'stiff': False}  # defaults in CC3D v3.7.7
###       SBMLoptions = {'relative': 1e-10, 'absolute': 1e-12, 'steps': 1000, 'stiff': True}
###       self.setSBMLGlobalOptions(SBMLoptions)    
        
        # add the PBPK in SBML model for APAP
        modelFile1 = 'Simulation/BIOMD0000000619.xml' # time unit=HOUR, substance=MOLE, volume=LITER, 
                                           # hasOnlySubstanceUnits="true", therfore unit is mass or moles
        initialConditions1={}
        initialConditions1['APAP_Dose_grams'] = pbpk_dose
        initialConditions1['BW'] = pbpk_bw
        initialConditions1['Fraction_unbound_plasma'] = pbpk_Fup
        initialConditions1['kGutabs'] = pbpk_kGutabs
        initialConditions1['Kkidney2plasma'] = pbpk_Kk2p
        initialConditions1['KRest2plasma'] = pbpk_Kr2p
        initialConditions1['Ratioblood2plasma'] = pbpk_Rb2p
        initialConditions1['Qgfr'] = pbpk_Qgfr
        initialConditions1['Fraction_unbound_plasma'] = pbpk_Fup
        initialConditions1['CLmetabolism'] = pbpk_CLmetabolism  ### turn off PBPK-SBML model's metabolism ODE
        initialConditions1['Kliver2plasma'] = pbpk_Kl2p  ###
        timeStep1=1./3600.*self.mcs2second # one second
        self.addFreeFloatingSBML(_modelFile=modelFile1,_modelName='PBPK_APAP',_stepSize=timeStep1,\
            _initialConditions=initialConditions1)    
 
        # add the PBPK in SBML model for APAP-Sulfate
        initialConditions2={}
        initialConditions2['APAP_Dose_grams'] = 0.
        initialConditions2['APAP_Dose'] = 0.
        initialConditions2['[AGutlumen]'] = 0.
        initialConditions2['AGutlumen'] = 0.
        initialConditions2['APAP_MW']= 231.2 # sulfate
        initialConditions2['BW'] = pbpk_bw
        initialConditions2['Fraction_unbound_plasma'] = pbpk_FupS
        initialConditions2['kGutabs'] = 0.
        initialConditions2['Kkidney2plasma'] = pbpk_Kk2pS
        initialConditions2['KRest2plasma'] = pbpk_Kr2pS
        initialConditions2['Ratioblood2plasma'] = pbpk_Rb2pS
        initialConditions2['Qgfr'] = pbpk_QgfrS
        initialConditions2['Fraction_unbound_plasma'] = pbpk_FupS
        initialConditions2['CLmetabolism'] = pbpk_CLmetabolism  # turn off PBPK-SBML model's metabolism ODE
        initialConditions2['Kliver2plasma'] = pbpk_Kl2pS  ###
        self.addFreeFloatingSBML(_modelFile=modelFile1,_modelName='PBPK_APAPS',_stepSize=timeStep1,\
            _initialConditions=initialConditions2)    

        # add the PBPK in SBML model for APAP-Glucuronide
        initialConditions3={}
        initialConditions3['APAP_Dose_grams'] = 0.
        initialConditions3['APAP_Dose'] = 0.
        initialConditions3['[AGutlumen]'] = 0.
        initialConditions3['AGutlumen'] = 0.
        initialConditions3['APAP_MW']= 327.3 # glucuronide
        initialConditions3['BW'] = pbpk_bw
        initialConditions3['Fraction_unbound_plasma'] = pbpk_FupG
        initialConditions3['kGutabs'] = 0.
        initialConditions3['Kkidney2plasma'] = pbpk_Kk2pG
        initialConditions3['KRest2plasma'] = pbpk_Kr2pG
        initialConditions3['Ratioblood2plasma'] = pbpk_Rb2pG
        initialConditions3['Qgfr'] = pbpk_QgfrG
        initialConditions3['Fraction_unbound_plasma'] = pbpk_FupG
        initialConditions3['CLmetabolism'] = pbpk_CLmetabolism  # turn off PBPK-SBML model's metabolism ODE
        initialConditions3['Kliver2plasma'] = pbpk_Kl2pG  ###
        self.addFreeFloatingSBML(_modelFile=modelFile1,_modelName='PBPK_APAPG',_stepSize=timeStep1,\
            _initialConditions=initialConditions3)    

        # add the subcellular metabolism model to the one hepatocyte
        modelFile2 = 'Simulation/BIOMD0000000624.xml' # time unit=SECONDS, substance=MILLIMOLE, volume=LITER, 
                                              # hasOnlySubstanceUnits="false", therfore unit is concentration
        initialConditions4={}
        initialConditions4['kGsh'] = sc_kGsh
        initialConditions4['Km_2E1_APAP'] = sc_Km_2E1_APAP
        initialConditions4['Km_PhaseIIEnzGlu_APAP'] = sc_Km_GLUC
        initialConditions4['Km_PhaseIIEnzSul_APAP'] = sc_Km_SULF
        initialConditions4['kNapqiGsh'] = sc_kNapqiGsh
        initialConditions4['Vmax_2E1_APAP'] = sc_Vmax_2E1_APAP
        initialConditions4['Vmax_PhaseIIEnzGlu_APAP'] = sc_Vmax_GLUC
        initialConditions4['Vmax_PhaseIIEnzSul_APAP'] = sc_Vmax_SULF
        initialConditions4['APAP']  = 0.
        timeStep2=1.*self.mcs2second # one second
        hep1=self.attemptFetchingCellById(1)
        self.addSBMLToCellIds(_modelFile=modelFile2,_modelName='METAB_APAP',_ids=[hep1.id],_stepSize=timeStep2,\
            _initialConditions=initialConditions4)       
        
        # returns dictionary of values for the single cell
        SUBCELLspeciesDict=self.getSBMLState(_modelName='METAB_APAP',_cell=hep1) 
        for key,value in SUBCELLspeciesDict.items(): 
            print "Start SUBCELL:\t",key,value
        PBPKspeciesDictAPAP =self.getSBMLState(_modelName='PBPK_APAP')
        PBPKspeciesDictAPAPG=self.getSBMLState(_modelName='PBPK_APAPG')
        PBPKspeciesDictAPAPS=self.getSBMLState(_modelName='PBPK_APAPS')
        for key,value in PBPKspeciesDictAPAP.items(): 
            print "Start PBPK A,G,S: ",key,value,PBPKspeciesDictAPAPG[key],PBPKspeciesDictAPAPS[key]
        print "\n\n"
            
        # setup the plot  TAKES 5X LONGER WITH THE PLOT WINDOWS EXPANDED!
        self.doPlot=True # True or False ###########################################
        if self.doPlot:
            self.pW1=self.addNewPlotWindow(_title='Plasma Concentrations',_xAxisTitle='Time (hours)', \
                _yAxisTitle='ug/ml',_xScaleType='linear',_yScaleType='linear')
            self.pW1.addPlot('APAP',     _style='Lines',_color='red',   _size=2)
            self.pW1.addPlot('APAPG',    _style='Lines',_color='green', _size=2)
            self.pW1.addPlot('APAPS',    _style='Lines',_color='yellow',_size=2)
            self.pW1.addPlot('APAP_exp', _style='Dots', _color='red',   _size=4)
            self.pW1.addPlot('APAPG_exp',_style='Dots', _color='green', _size=4)
            self.pW1.addPlot('APAPS_exp',_style='Dots', _color='yellow',_size=4)
            # plot the t=0 points
            self.pW1.addDataPoint("APAP", 0,0); self.pW1.addDataPoint("APAPS",0,0); self.pW1.addDataPoint("APAPG",0,0)
            # Plot the human ADME data of Critchley et al, data is ug/ml in plasma
            # Data:time(hr),APAP,APAPS,APAPG
            self.ADME=( (0.0, 0.0, 0.0,  0.0),  \
                        (0.3, 8.94,0.76, 0.33), ( 0.5,14.67,1.97, 1.75), (0.75,16.77,2.96,3.65), \
                        (1.0,16.63,3.41, 4.71), ( 1.5,14.89,4.06, 7.74), (2.0, 13.98,4.24,9.32), \
                        (3.0,10.85,4.24,10.78), ( 4.0, 8.14,3.71,10.75), (5.0,  6.12,3.20,9.48), \
                        (6.0, 4.45,2.68, 8.13), ( 7.0, 3.13,2.15, 6.78), (8.0,  2.72,1.66,5.18), \
                       (12.0, 1.20,0.72, 2.05), (24.0, 0.38,0.12, 0.27) )
            self.ADMEtimeList =[x[0] for x in self.ADME]
            self.ADMEout = "" # will be a space delimited text 
            self.ADMEout = " MCS  Time(hr)  " \
                + "APAP_exp_ug/ml  APAP_calc_ug/ml  APAP_abs_err_ug/ml  "    \
                + "APAPS_exp_ug/ml  APAPS_calc_ug/ml  APAPS_abs_err_ug/ml  " \
                + "APAPG_exp_ug/ml  APAPG_calc_ug/ml  APAPG_abs_err_ug/ml  " \
                + "total_err_this_step  cumulative_total_err  "              \
                + "RMSeTotal  accumChiSquaredTotal\n"

            self.pW2=self.addNewPlotWindow(_title='Total Recovered APAP Equivalent (in Tubules Only)',\
                _xAxisTitle='Time (hours)',_yAxisTitle='Dose APAP Equiv',_xScaleType='linear',_yScaleType='linear')
            self.pW2.addPlot('Total', _style='Lines',_color='Red',   _size=1)
            self.pW2.addPlot('Theory',_style='Lines',_color='Green', _size=2)

            self.pW3=self.addNewPlotWindow(_title='In Hep1', \
                _xAxisTitle='Time (hours)',_yAxisTitle='milliMolarity', \
                _xScaleType='linear',_yScaleType='linear')
            self.pW3.addPlot('APAP',    _style='Lines', _color='Red',   _size=2)
            self.pW3.addPlot('APAPG',   _style='Lines', _color='Green', _size=2)
            self.pW3.addPlot('APAPS',   _style='Lines', _color='Yellow',_size=2)
            self.pW3.addPlot('NAPQIGSH',_style='Dots',  _color='Blue',  _size=2)
            
            # some accumulators for the statistics
            self.cumTotErr = 0.
            self.accumRMSeTotal = 0.
            self.accumChiSquaredTotal = 0.
            self.matchedPtCount = 0.
        
    ##########################################################################################
    def step(self,mcs):
        if mcs == 1:
            # plot the in vivo ADME data
            for t,A,S,G in self.ADME:  # this doesn't work in the init!
                if t <= totalSimTimeHr:  # only plot the data needed given the simulation time
                    self.pW1.addDataPoint("APAP_exp", t,A)
                    self.pW1.addDataPoint("APAPS_exp",t,S)
                    self.pW1.addDataPoint("APAPG_exp",t,G)
            # plot the Theoretical limit points
            self.pW2.addDataPoint("Theory", 0.,1.)
            self.pW2.addDataPoint("Theory", 8.,1.)

        if mcs % 1 == 0:        
            # time step all the SBML models
            self.timestepSBML()
            hep1  =self.attemptFetchingCellById(1)  # the first hep
            blood1=self.attemptFetchingCellById(2)  # the first blood "cell"
            
            # get the values from the three PBPK models for APAP, APAPS and APAPG
            # Units are MOLES, need to convert to milliMolar using the liver's PBPK volume
            PBPKspeciesDictAPAP =self.getSBMLState(_modelName='PBPK_APAP')
            PBPKspeciesDictAPAPS=self.getSBMLState(_modelName='PBPK_APAPS')
            PBPKspeciesDictAPAPG=self.getSBMLState(_modelName='PBPK_APAPG')
            # get the subcellular values
            SUBCELLspeciesDict  =self.getSBMLState(_modelName='METAB_APAP',_cell=hep1) # dict of values for the one cell

            # update the CC3D blood compartment/cell from the PBPK model
            # the values returned from the SBML PBPK are in Moles!
            # convert to millimole/liter       vvvvvvvvvvvvvvvvvvvvvv 1/25 vvvvvvvvvvvvvvvvvvvv
            # APAP concentration (millimole/Liter) in liver compartment
            blood1.dict['APAP']  = PBPKspeciesDictAPAP['CLiver'] *1000./PBPKspeciesDictAPAP['VLiver'] 
            blood1.dict['APAPS'] = PBPKspeciesDictAPAPS['CLiver']*1000./PBPKspeciesDictAPAPS['VLiver']
            blood1.dict['APAPG'] = PBPKspeciesDictAPAPG['CLiver']*1000./PBPKspeciesDictAPAPG['VLiver']

            # update the Hep from the subcellular model
            # SBML is mM concentration
            hep1.dict['APAP'] =SUBCELLspeciesDict['APAP']  # SBML is mM concentration
            hep1.dict['APAPG']=SUBCELLspeciesDict['APAPconj_Glu']
            hep1.dict['APAPS']=SUBCELLspeciesDict['APAPconj_Sul']
            hep1.dict['NAPQIGSH']=SUBCELLspeciesDict['NAPQIGSH']  
 
            ###4####################
            # NOTE that in the Liver Paper I the hep <--> blood transfers fluxes are scaled by the common surface area
            # but not doing that here
            #
            bloodHepSurfaceArea = 20 # um^2? pixel   ###4  this is the max surface area (actually length) in the Liver I model
            bloodHepSurfaceArea = 20*10*2. # um^2? pixel   ###5  this is the max surface area (actually length) in the Liver I model
            bloodHepSurfaceArea = 1.0 # 
            #
            # in the CC3D model the time unit was seconds and amounts in mMol/liter  ################
            #
            #  hepatocyte volume =8x10^-15L  (20um)^3
            #  sinusoid segment volume is 1.3x10^-15L (14% of hep + sinuoid)
            #  ratio = 8/1.3 = 6.15, which is the volume ratio of heps to sinusoid
            hepOverBloodVolumeRatio = hepatoVolume/bloodVolume
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            hepOverBloodVolumeRatio = 1.
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ############################## START iterate CC3D ODEs #####################################
            for iii in range(1,cc3d_iterODE+1):  # extra steps for the CC3D ODE transfer calculations

                # Calculate passive diffusion, only blood <=> hep
                # cc3d_k_PD_S2H # using the Serum <--> Hep constants, Fup (need Hc? and RB2P?)
                hep2bloodFlux = cc3d_k_PD_H2S*hep1.dict['APAP']               \
                        * bloodHepSurfaceArea *self.mcs2second /cc3d_iterODE ###4
                blood2hepFlux = cc3d_k_PD_S2H*blood1.dict['APAP'] * pbpk_Fup/pbpk_Rb2p  \
                    * bloodHepSurfaceArea *self.mcs2second  /cc3d_iterODE ###4
                hep1.dict['APAP']   -= hep2bloodFlux 
                hep1.dict['APAP']   += blood2hepFlux /hepOverBloodVolumeRatio
                blood1.dict['APAP'] += hep2bloodFlux /hepOverBloodVolumeRatio
                blood1.dict['APAP'] -= blood2hepFlux
               
                # Calculate APAP active import blood -> hep
                # cc3d_Km_AT_APAP, cc3d_Vmax_AT_APAP,  this uses Fup  (need Hc? and RB2P?)
                # NOTE that the Km and Vmax are in mM
                APAPactFlux = cc3d_Vmax_AT_APAP*blood1.dict['APAP']*pbpk_Fup/pbpk_Rb2p \
                    / (cc3d_Km_AT_APAP + blood1.dict['APAP']*pbpk_Fup/pbpk_Rb2p)       \
                    * bloodHepSurfaceArea *self.mcs2second  /cc3d_iterODE  ###4
                hep1.dict['APAP']   += APAPactFlux /hepOverBloodVolumeRatio
                blood1.dict['APAP'] -= APAPactFlux 
    
                # Calculate the active export of APAPG and APAPS from the Hep -> Blood
                METAB_G_actFlux = cc3d_k_AT_APAPG*hep1.dict['APAPG'] * \
                    bloodHepSurfaceArea *self.mcs2second /cc3d_iterODE ###4
                METAB_S_actFlux = cc3d_k_AT_APAPS*hep1.dict['APAPS'] * \
                    bloodHepSurfaceArea *self.mcs2second /cc3d_iterODE ###4
                blood1.dict['APAPG'] += METAB_G_actFlux *hepOverBloodVolumeRatio
                blood1.dict['APAPS'] += METAB_S_actFlux *hepOverBloodVolumeRatio
                hep1.dict['APAPG']   -= METAB_G_actFlux
                hep1.dict['APAPS']   -= METAB_S_actFlux

#                 print "mcs,iii,hep2bloodFlux,blood2hepFlux,APAPactFlux,METAB_G_actFlux,METAB_S_actFlux:\n", \
#                     mcs,iii,hep2bloodFlux,blood2hepFlux,APAPactFlux,METAB_G_actFlux,METAB_S_actFlux
#                 print "blood1.dict['APAP'],hep1.dict['APAP']:\n",blood1.dict['APAP'],hep1.dict['APAP']

            ############################## END iterate CC3D ODEs #####################################

            # update the subcellular model for the changes in CC3D APAP/APAPG/APAPS  
            self.setSBMLValue(_modelName='METAB_APAP',_valueName='APAP',        _value=hep1.dict['APAP'] ,_cell=hep1)
            self.setSBMLValue(_modelName='METAB_APAP',_valueName='APAPconj_Glu',_value=hep1.dict['APAPG'],_cell=hep1)
            self.setSBMLValue(_modelName='METAB_APAP',_valueName='APAPconj_Sul',_value=hep1.dict['APAPS'],_cell=hep1)

            # put the value from the CC3D back into PBPKs
            # APAP amount (moles) in liver compartment 
            self.setSBMLValue(_modelName='PBPK_APAP', _valueName='CLiver',_value=blood1.dict['APAP'] /1000.*PBPKspeciesDictAPAP['VLiver'])
            self.setSBMLValue(_modelName='PBPK_APAPG',_valueName='CLiver',_value=blood1.dict['APAPG']/1000.*PBPKspeciesDictAPAPG['VLiver'])
            self.setSBMLValue(_modelName='PBPK_APAPS',_valueName='CLiver',_value=blood1.dict['APAPS']/1000.*PBPKspeciesDictAPAPS['VLiver'])
 
        if mcs % 60 == 0:
            print "\n %6i AGS  %7.4f %7.4f %7.4f uM" % (mcs,hep1.dict['APAP']*1.e3,hep1.dict['APAPG']*1.e3,hep1.dict['APAPS']*1.e3)
            print "          %9.6f %9.6f\n" % (PBPKspeciesDictAPAP['CLiver']/PBPKspeciesDictAPAP['VLiver'],hep1.dict['APAP'])
            # update the cell display fields
            for cell in self.cellList:   
                self.scalarCL_APAP[cell] =cell.dict['APAP']
                self.scalarCL_APAPG[cell]=cell.dict['APAPG']
                self.scalarCL_APAPS[cell]=cell.dict['APAPS']
                
            # write to the log file # # # # # # # # # # # # # # # # # # # # # # # #
            outline  = str(mcs)
            outline += "\t"+str(mcs/60./60.*self.mcs2second)    # hours
            outline += "\t"+str(hep1.dict['APAP'])+"\t"+str(hep1.dict['APAPG'])+"\t"+str(hep1.dict['APAPS'])
            outline += "\t"+str(blood1.dict['APAP'])+"\t"+str(blood1.dict['APAPG'])+"\t"+str(blood1.dict['APAPS'])
            outline += "\t"+str(PBPKspeciesDictAPAP['CLiver'])
            outline += "\t"+str(PBPKspeciesDictAPAPG['CLiver'])
            outline += "\t"+str(PBPKspeciesDictAPAPS['CLiver'])
            outline += "\t"+str(PBPKspeciesDictAPAP['CVen'])
            outline += "\t"+str(PBPKspeciesDictAPAPG['CVen'])
            outline += "\t"+str(PBPKspeciesDictAPAPS['CVen'])
            outline += "\t"+str(PBPKspeciesDictAPAP['CTubules'])
            outline += "\t"+str(PBPKspeciesDictAPAPG['CTubules'])
            outline += "\t"+str(PBPKspeciesDictAPAPS['CTubules'])
            outline += "\t"+str(SUBCELLspeciesDict['APAP'])
            outline += "\t"+str(SUBCELLspeciesDict['APAPconj_Glu'])
            outline += "\t"+str(SUBCELLspeciesDict['APAPconj_Sul'])
            outline += "\t"+str(SUBCELLspeciesDict['NAPQI'])
            outline += "\t"+str(SUBCELLspeciesDict['NAPQIGSH'])
            outline += "\t"+str(SUBCELLspeciesDict['GSH'])+"\n"
            self.fileHandleLog.write(outline)

        if (mcs in (0,1,2,3,4)) or (mcs % 3600 == 0):
            print "\n\n"
            for key,value in SUBCELLspeciesDict.items(): 
                print "MCS=",mcs,"  SUBCELL:\t",key,value
            for key,value in PBPKspeciesDictAPAP.items(): 
                print "MCS=",mcs,"   PBPK A,G,S: ",key,value,PBPKspeciesDictAPAPG[key],PBPKspeciesDictAPAPS[key]
            print "\n\n"                   
            
        # update the plots, note that the plot y-axis is ug/ml (mg/L) but SBML-PBPK data is moles and it is 
        # mMolar in the cell Dict. 
        # X-axis in hours
        if self.doPlot & (mcs % 300 == 0) :
            aTime = mcs/3600.*self.mcs2second
            #                                       vvvvvvvvvvvvvvvvvvvvvv 1/25 vvvvvvvvvvvvvvvvv
            self.APAPug_ml  = PBPKspeciesDictAPAP['CVen'] *151.2*1000./PBPKspeciesDictAPAP['VVen']
            self.APAPSug_ml = PBPKspeciesDictAPAPS['CVen']*231.2*1000./PBPKspeciesDictAPAPS['VVen']
            self.APAPGug_ml = PBPKspeciesDictAPAPG['CVen']*327.3*1000./PBPKspeciesDictAPAPG['VVen']

            self.pW1.addDataPoint("APAP", aTime,self.APAPug_ml ) 
            self.pW1.addDataPoint("APAPG",aTime,self.APAPGug_ml) 
            self.pW1.addDataPoint("APAPS",aTime,self.APAPSug_ml) 

            totAPAPequivTubMoles = PBPKspeciesDictAPAP['CTubules'] + \
                PBPKspeciesDictAPAPG['CTubules'] + PBPKspeciesDictAPAPS['CTubules']
            fractTotAPAP = totAPAPequivTubMoles/(pbpk_dose/151.2)
            print "         CTubules A,G,S, fractOfTotal:",PBPKspeciesDictAPAP['CTubules'], \
                PBPKspeciesDictAPAPG['CTubules'],PBPKspeciesDictAPAPS['CTubules'],fractTotAPAP
            self.pW2.addDataPoint("Total",aTime,fractTotAPAP)
            
            self.pW3.addDataPoint('APAP',    aTime,hep1.dict['APAP'])
            self.pW3.addDataPoint('APAPG',   aTime,hep1.dict['APAPG'])
            self.pW3.addDataPoint('APAPS',   aTime,hep1.dict['APAPS'])
            self.pW3.addDataPoint('NAPQIGSH',aTime,hep1.dict['NAPQIGSH'])
                        
        # Calculate errors
        # using the in vivo plot data, convert time, conc etc.
        timeHr = mcs/60./60.*self.mcs2second
        if timeHr in self.ADMEtimeList:  # the current mcs translates into one of the times in the ADME data list
            print "\n\n"
            i=self.ADMEtimeList.index(timeHr)
            self.matchedPtCount +=3  #  +3 for APAP, APAPG, APAPS
            dAPAP =abs(self.ADME[i][1]-self.APAPug_ml)
            dAPAPS=abs(self.ADME[i][2]-self.APAPSug_ml)
            dAPAPG=abs(self.ADME[i][3]-self.APAPGug_ml)
            dTotal=dAPAP+dAPAPG+dAPAPS # in ug/ml
            self.cumTotErr += dTotal
            self.accumRMSeTotal += dAPAP**2+dAPAPS**2+dAPAPG**2
            self.RMSeTotal = math.sqrt(self.accumRMSeTotal/self.matchedPtCount)
            if self.ADME[i][1] <> 0. and self.ADME[i][2] <> 0. and self.ADME[i][3] <> 0. :
                chiSquared =  dAPAP**2/self.ADME[i][1] + dAPAPS**2/self.ADME[i][2] + dAPAPG**2/self.ADME[i][3]
                self.accumChiSquaredTotal += chiSquared

            self.ADMEout += '%7i  %6.2f   '           % (mcs,timeHr)
            self.ADMEout += '%11.2f  %7.3f  %7.3f   ' % (self.ADME[i][1],self.APAPug_ml, dAPAP) 
            self.ADMEout += '%11.2f  %7.3f  %7.3f   ' % (self.ADME[i][2],self.APAPSug_ml,dAPAPS) 
            self.ADMEout += '%11.2f  %7.3f  %7.3f   ' % (self.ADME[i][3],self.APAPGug_ml,dAPAPG) 
            self.ADMEout += '%11.3f  %7.3f   '        % (dTotal,self.cumTotErr)
            self.ADMEout += '%11.4f    '              % (self.RMSeTotal)
            self.ADMEout += '%11.4f\n'                % (self.accumChiSquaredTotal)
            print self.ADMEout

    def finish(self):
        # Finish Function gets called after the last MCS
        # save the data for the error calc
        try:                
            fileHandle0,fullFileName0=self.openFileInSimulationOutputDirectory("ADME_error_data.txt","w")
        except IOError:
            print "Could not open file for writing ADME error data. "                
            return
        #f = open('workfile', 'w')
        fileHandle0.write(self.ADMEout)
        fileHandle0.close()
        
        self.fileHandleLog.close()
        
        # save the calculated ADME data
        if self.doPlot:
            try:                
                fileHandle1,fullFileName1=self.openFileInSimulationOutputDirectory("ADME_plot.png","w")
            except IOError:
                print "Could not open file for writing ADME image. "                
                return
            try:                
                fileHandle2,fullFileName2=self.openFileInSimulationOutputDirectory("ADME_plot.txt","w")
            except IOError:
                print "Could not open file for writing ADME data. "                
                return
            self.pW1.savePlotAsPNG(fullFileName1,1000,1000) # image size default is 400 x 400       
            self.pW1.savePlotAsData(fullFileName2) 
            
        #####################################################################
        # write a summary for parameter scanning, list all the parameters in the parameter file
        try:                
            aFile,aFilefullFileName=self.openFileInSimulationOutputDirectory("scan_summary.txt","w")
        except IOError:
            print "Could not open file for writing summary data. "                
            return
        #aFilefullFileName ="c:/scan_summary.txt"   # use to put all output into a single file
        #aFile = open(aFilefullFileName,'a+')
        
        outline1 ="";   outline2=""
#        for aVar in dir(p):   # dir returns a sorted list
        for aVar in dir():   # dir returns a sorted list
            if not aVar.startswith('__'):  # skip the special variaiables and objects like __init__
                #print aVar,"\t",__dict__.get(aVar),"\n"
                outline1 += aVar+"\t"  # headers line, tab delimited
                outline2 += str(__dict__.get(aVar))+"\t"  # the value of aVar, values line, tab delimited
       
        outline1 +="self.cumTotErr"+"\t"
        outline2 +=str(self.cumTotErr)+"\t"
        outline1 +="self.RMSeTotal"+"\t"
        outline2 +=str(self.RMSeTotal)+"\t"
        outline1 +="self.accumChiSquaredTotal"+"\t"
        outline2 +=str(self.accumChiSquaredTotal)+"\t"
        outline1 +="FileName"+"\t"
        outline2 +=aFilefullFileName+"\t"
        outline1 += "\n"
        outline2 += "\n"
       
        aFile.write(outline1)
        aFile.write(outline2)
        aFile.close()
