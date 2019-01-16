# total simulation time 
totalSimTimeHr = 12 # 8 hours default, in vitro data goes out to 24 hours

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
pbpk_CLmetabolism=0.0    # 1/hr

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

cc3d_k_AT_APAPG	= 0.00005 # 1/s  0.000450  1/30
cc3d_k_AT_APAPS	= 0.00040 # 1/s  0.001900

cc3d_Km_AT_APAP	 = 0.00200	 # mmol/L  0.010
cc3d_Vmax_AT_APAP= 0.00060	 # mmol/L/s  0.010  0.005 0.00010 0.00030 #######1/26/18#######

# cc3d_k_PD_H2S	= 0.000050       # 1/s   0.00100
# cc3d_k_PD_S2H	= cc3d_k_PD_H2S  # 1/s  0.00100
