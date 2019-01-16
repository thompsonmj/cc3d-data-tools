# total simulation time 
totalSimTimeHr = 24 # 8 hours default, in vitro data goes out to 24 hours

# number of iterations per mcs on the CC3D's model ODEs
cc3d_iterODE = 1  # best = 5, Default = 1. If 5, then CC3D will divide the transfer ODE calc's into 5 substeps

# conversion factor for CC3D mcs to simulated time
# note that if this value is changed the total number of mcs needs to be changed as well
cc3d_mcs2second = 1 # Best = 0.5, Default = 1.0 simulated second per mcs (second/mcs), if 0.5 then 2 mcs per simulated second


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
hepatoVolume = 8.0e-12 # Liter, 20.**3*1e-18*1000. the compartment is one cell, 20um^3
bloodVolume  = 1.3e-12 # Liter, this is the blood volume in a hep long sinusoid,
                       # based on ~15% of liver volume being sinusoid			
 
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
pbpk_dose	= 1.400	 # g       1.4
pbpk_hemat	= 0.450	 # unitless
pbpk_kGutabs	= 1.400	 # 1/hr		
pbpk_CLmetabolism=0.0    # 1/hr

cc3d_Km_AT_APAP	= 0.010	 # mmol/L  0.010
cc3d_Vmax_AT_APAP=0.0003	 # mmol/L/s  0.010   xxxxxxx
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
sc_Vmax_2E1_APAP= 0.000020 # mmol/L/s   
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
 
########################################################################
########################################################################

pbpk_kGutabs	= 1.4	 # 1/hr	                        1.5,1.7,1.9

cc3d_k_AT_APAPG	= 0.0006 # 1/s  0.000450               5e-05,0.00025,0.00045
cc3d_k_AT_APAPS	= 0.0026 # 1/s  0.001900               0.0004,0.00115,0.0019

cc3d_Km_AT_APAP	 = 0.01500 # mmol/L    0.010            0.0001,0.00505,0.01
cc3d_Vmax_AT_APAP= 0.00040 # mmol/L/s  0.010            0.010,0.005,0.00010,0.00030

sc_Km_GLUC	= 1.500	   # mmol/L
sc_Vmax_GLUC	= 0.00120  # mmol/L/s                   0.0003,0.0006,0.0009,0.0012

sc_Km_SULF	= 0.200	   # mmol/L
sc_Vmax_SULF	= 0.000150 # mmol/L/s                   5e-05,0.0001,0.00015,0.0002

cc3d_k_PD_H2S	= 0.00150       # 1/s   0.00100        2.5e-05,5e-05,0.00010
cc3d_k_PD_S2H	= cc3d_k_PD_H2S  # 1/s  0.00100

######################## ChiSq=107 ################################################

pbpk_kGutabs	= 1.0	 # 1/hr	                        1.5,1.7,1.9

cc3d_k_AT_APAPG	= 0.001 # 1/s  0.000450               5e-05,0.00025,0.00045
cc3d_k_AT_APAPS	= 0.005 # 1/s  0.001900               0.0004,0.00115,0.0019

cc3d_Km_AT_APAP	 = 0.0200 # mmol/L    0.010            0.0001,0.00505,0.01
cc3d_Vmax_AT_APAP= 0.00050 # mmol/L/s  0.010            0.010,0.005,0.00010,0.00030

sc_Km_GLUC	= 1.630	   # mmol/L
sc_Vmax_GLUC	= 0.00120  # mmol/L/s                   0.0003,0.0006,0.0009,0.0012

sc_Km_SULF	= 0.100	   # mmol/L
sc_Vmax_SULF	= 0.00020 # mmol/L/s                   5e-05,0.0001,0.00015,0.0002

cc3d_k_PD_H2S	= 0.00150       # 1/s   0.00100        2.5e-05,5e-05,0.00010
cc3d_k_PD_S2H	= cc3d_k_PD_H2S  # 1/s  0.00100

######################## ChiSq=141 ################################################
pbpk_kGutabs	= 1.4	 # 1/hr	                        1.5,1.7,1.9

cc3d_k_AT_APAPG	= 0.0006 # 1/s  0.000450               5e-05,0.00025,0.00045
cc3d_k_AT_APAPS	= 0.0026 # 1/s  0.001900               0.0004,0.00115,0.0019

cc3d_Km_AT_APAP	 = 0.01500 # mmol/L    0.010            0.0001,0.00505,0.01
cc3d_Vmax_AT_APAP= 0.0010 # mmol/L/s  0.010  ***        0.010,0.005,0.00010,0.00030

sc_Km_GLUC	= 1.500	   # mmol/L
sc_Vmax_GLUC	= 0.00120  # mmol/L/s                   0.0003,0.0006,0.0009,0.0012

sc_Km_SULF	= 0.200	   # mmol/L
sc_Vmax_SULF	= 0.000150 # mmol/L/s                   5e-05,0.0001,0.00015,0.0002

cc3d_k_PD_H2S	= 0.00150       # 1/s   0.00100        2.5e-05,5e-05,0.00010
cc3d_k_PD_S2H	= cc3d_k_PD_H2S  # 1/s  0.00100

######################## ChiSq=115 ################################################
pbpk_kGutabs	= 1.4	 # 1/hr	                        1.5,1.7,1.9

cc3d_k_AT_APAPG	= 0.0009 # 1/s  0.000450      ***         5e-05,0.00025,0.00045
cc3d_k_AT_APAPS	= 0.005 # 1/s  0.001900       ***        0.0004,0.00115,0.0019

cc3d_Km_AT_APAP	 = 0.01500 # mmol/L    0.010            0.0001,0.00505,0.01
cc3d_Vmax_AT_APAP= 0.0008 # mmol/L/s  0.010  ***        0.010,0.005,0.00010,0.00030

sc_Km_GLUC	= 1.500	   # mmol/L
sc_Vmax_GLUC	= 0.00120  # mmol/L/s                   0.0003,0.0006,0.0009,0.0012

sc_Km_SULF	= 0.200	   # mmol/L
sc_Vmax_SULF	= 0.000150 # mmol/L/s                   5e-05,0.0001,0.00015,0.0002

cc3d_k_PD_H2S	= 0.00150       # 1/s   0.00100        2.5e-05,5e-05,0.00010
cc3d_k_PD_S2H	= cc3d_k_PD_H2S  # 1/s  0.00100

# # ######################## ChiSq=44 ################################################
# # pbpk_QgfrG	= 12.	 # L/hr   7.86  12 15 18 22##############
# # pbpk_QgfrS	= 12.	 # L/hr   9.96  15 20 22  ##############
# # ######################## ChiSq=28.7 ##############################################
# # pbpk_QgfrG	= 14.	 # L/hr   7.86  12 15 18 22##############
# # pbpk_QgfrS	= 14.	 # L/hr   9.96  15 20 22  ##############
# # ######################## ChiSq=15.6 ##############################################
# # pbpk_QgfrG	= 18.	 # L/hr   7.86  12 15 18 22##############
# # pbpk_QgfrS	= 18.	 # L/hr   9.96  15 20 22  ##############
# # ######################## ChiSq=12.9 ##############################################
# # pbpk_QgfrG	= 22.	 # L/hr   7.86  12 15 18 22##############
######################## ChiSq=10.1 ##############################################
pbpk_QgfrG	= 20.	 # L/hr   7.86  12 15 18 22##############
pbpk_QgfrS	= 18.	 # L/hr   9.96  15 20 22  ##############
cc3d_k_AT_APAPG	= 0.00045 #1/s  0.0006

# ######################## ChiSq=11.03 ################################################
# pbpk_kGutabs	= 1.2	 # 1/hr	                        1.5,1.7,1.9
# ######################## ChiSq=10.3  ################################################
# pbpk_kGutabs	= 1.6	 # 1/hr	                        1.5,1.7,1.9

# ######################## ChiSq=32 ################################################
# cc3d_Km_AT_APAP	 = 0.0100 # mmol/L    0.015            0.0001,0.00505,0.01
# cc3d_Vmax_AT_APAP= 0.0010 # mmol/L/s  0.0008  ***        0.010,0.005,0.00010,0.00030
# ######################## ChiSq=4.7 ################################################
cc3d_Km_AT_APAP	 = 0.0200 # mmol/L    0.015            0.0001,0.00505,0.01
cc3d_Vmax_AT_APAP= 0.0006 # mmol/L/s  0.0008  ***        0.010,0.005,0.00010,0.00030

# ######################## ChiSq=5.85 ################################################
# sc_Km_GLUC	= 1.600	   # mmol/L    1.50
# ######################## ChiSq=3.76 ################################################
# sc_Km_GLUC	= 1.400	   # mmol/L
# ######################## ChiSq=3.18 ################################################
# sc_Km_GLUC	= 1.300	   # mmol/L
# ######################## ChiSq=3.00 ################################################
sc_Km_GLUC	= 1.200	   # mmol/L

# ######################## ChiSq=4.46 ################################################
# sc_Vmax_GLUC= 0.00100	# mmol/L/s	0.0009    was 0.0012      
# ######################## ChiSq=4.18 ################################################
#sc_Vmax_GLUC= 0.00140	# mmol/L/s	0.0009    was 0.0012      

# ######################## ChiSq=2.91 ################################################
# cc3d_Km_AT_APAP = 0.0250 # mmol/L    0.010  was 0.02          
# ######################## ChiSq=3.87 ################################################
# cc3d_Km_AT_APAP	 = 0.0300 # mmol/L    0.010  was 0.02          


#=================================================================================
######################## ChiSq=2.92 ################################################
# pbpk_kGutabs	= 1.4	 # 1/hr	    1.4                  

# cc3d_k_AT_APAPG	= 0.00045 # 1/s  0.000450              
# cc3d_k_AT_APAPS	= 0.005 # 1/s  0.001900            

# cc3d_Km_AT_APAP	 = 0.0250 # mmol/L     0.010           
# cc3d_Vmax_AT_APAP= 0.0006 # mmol/L/s   0.010          

# sc_Km_GLUC	= 1.200	   # mmol/L    1.5
# sc_Vmax_GLUC	= 0.00120  # mmol/L/s  0.0009         

# sc_Km_SULF	= 0.200	   # mmol/L    0.3
# sc_Vmax_SULF	= 0.000150 # mmol/L/s  0.00015              

# cc3d_k_PD_H2S	= 0.00150       # 1/s   0.00100        
# cc3d_k_PD_S2H	= cc3d_k_PD_H2S  # 1/s  0.00100

# pbpk_QgfrG	= 20.	 # L/hr   7.86  
# pbpk_QgfrS	= 18.	 # L/hr   9.96  

#========================= Pscan 7 ===========================================
######################## ChiSq=2.17 ##########################################
pbpk_kGutabs	= 1.4	 # 1/hr	    1.4    

cc3d_k_AT_APAPG	= 0.00046 # 1/s  0.000450              
cc3d_k_AT_APAPS	= 0.00495 # 1/s  0.001900            

cc3d_Km_AT_APAP	 = 0.0250 # mmol/L     0.010           
cc3d_Vmax_AT_APAP= 0.00046 # mmol/L/s   0.010          

sc_Km_GLUC	= 1.200	   # mmol/L    1.5
sc_Vmax_GLUC	= 0.00162  # mmol/L/s  0.0009         

sc_Km_SULF	= 0.200	   # mmol/L    0.3
sc_Vmax_SULF	= 0.000215 # mmol/L/s  0.00015              

cc3d_k_PD_H2S	= 0.00125       # 1/s   0.00100        
cc3d_k_PD_S2H	= cc3d_k_PD_H2S  # 1/s  0.00100

pbpk_QgfrG	= 20.0	 # L/hr   7.86  
pbpk_QgfrS	= 21.8	 # L/hr   9.96  

# cc3d_Vmax_AT_APAP =	0.000460
# cc3d_k_AT_APAPG	=	0.000495
# cc3d_k_AT_APAPS	=	0.00450
# cc3d_k_PD_H2S	=	0.00125
# pbpk_QgfrG	=	20.0
# pbpk_QgfrS	=	21.8
# sc_Vmax_GLUC	=	0.00162
# sc_Vmax_SULF	=	0.000215

##### from particle swarm ending 2/26/18, ChiSq=10.07 ####################
pbpk_kGutabs	= 1.064
cc3d_k_AT_APAPG	= 0.0033
cc3d_k_AT_APAPS	= 0.0102
cc3d_Km_AT_APAP	= 0.239
cc3d_Vmax_AT_APAP= 0.0033
sc_Km_GLUC	= 7.5368
sc_Vmax_GLUC	= 0.0122
sc_Km_SULF	= 1.4067
sc_Vmax_SULF	= 0.0018
cc3d_k_PD_S2H	= 0.0036
cc3d_k_PD_H2S	= 0.0036     # added
pbpk_QgfrG	= 14.03
pbpk_QgfrS	= 38.66
BestE	        = 10.07

##### from particle swarm ending 2/28/18, ChiSq=2.01 iter=191 swrm=1 partcl=14 #############
pbpk_kGutabs	= 1.75
cc3d_k_AT_APAPG	= 1.17e-03
cc3d_k_AT_APAPS	= 3.17e-02
cc3d_Km_AT_APAP	= 1.49e-01
cc3d_Vmax_AT_APAP= 2.99e-03
sc_Km_GLUC	= 11
sc_Vmax_GLUC	= 1.83e-03
sc_Km_SULF	= 0.603
sc_Vmax_SULF	= 1.99e-04
cc3d_k_PD_S2H	= 8.08e-04
cc3d_k_PD_H2S	= 8.08e-04     # added
pbpk_QgfrG	= 11.2
pbpk_QgfrS	= 38.3
BestE	        = 2.01
#1.75e+00  1.17e-03  3.17e-02  1.49e-01  2.99e-03  1.10e+01  1.83e-03  6.03e-01  1.99e-04  8.08e-04  1.12e+01  3.83e+01  2.0065e+00    1      1       1   

##### from particle swarm ending 2/28/18, ChiSq=1.998 iter=198 swrm=1 partcl=17 #############
pbpk_kGutabs	= 1.74
cc3d_k_AT_APAPG	= 1.27e-03
cc3d_k_AT_APAPS	= 3.17e-02
cc3d_Km_AT_APAP	= 1.53e-01
cc3d_Vmax_AT_APAP= 3.00e-03
sc_Km_GLUC	= 10.9
sc_Vmax_GLUC	= 1.83e-03
sc_Km_SULF	= 0.603
sc_Vmax_SULF	= 2.03e-04
cc3d_k_PD_S2H	= 8.06e-04
cc3d_k_PD_H2S	= 8.06e-04     # added
pbpk_QgfrG	= 11.2
pbpk_QgfrS	= 38.2
BestE	        = 2.01
#198 1 17 1.74e+00 1.27e-03 3.17e-02 1.53e-01 3.00e-03 1.09e+01 1.83e-03 6.03e-01 2.03e-04 8.06e-04 1.12e+01 3.82e+01 1.9975 1 1 1   
