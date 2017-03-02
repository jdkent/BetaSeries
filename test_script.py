import os
outDir=os.path.join(os.getcwd(),'BetaSeries_testout')
func=os.path.join(os.getcwd(),'pre_sub10_FLANKER_RPI.nii.gz')
t1brain=os.path.join(os.getcwd(),'pre_sub10_MPRAGE-1_RPI_ss.nii.gz')
t1head=os.path.join(os.getcwd(),'pre_sub10_MPRAGE-1_RPI.nii.gz')
import BetaSeries_functions as BS
ICA_inputs=BS.PreICA(outDir=outDir,func=func,T1_brain=t1brain,T1_head=t1head)
fslDir=os.path.join(os.environ["FSLDIR"],'bin','')
ICA_AROMAoutDir=os.path.join(outDir,"ICA_AROMA")
melDir=os.path.join(ICA_AROMAoutDir,"melodic.ica")
TR=2
dim=0
denType='nonaggr'
melIC = os.path.join(ICA_AROMAoutDir,'melodic_IC_thr.nii.gz')
melIC_MNI =  os.path.join(ICA_AROMAoutDir,'melodic_IC_thr_MNI2mm.nii.gz')
if not os.path.isdir(ICA_AROMAoutDir):
	os.makedirs(ICA_AROMAoutDir)
BS.runICA(fslDir, ICA_inputs.func, ICA_AROMAoutDir, melDir, ICA_inputs.func_mask, dim, TR
BS.register2MNI(fslDir, melIC, melIC_MNI, ICA_inputs.functoT1_transform, ICA_inputs.T1toMNI_transform)
edgeFract, csfFract = BS.feature_spatial(fslDir, ICA_AROMAoutDir, '/home/james/dev/BetaSeries/', melIC_MNI)
melmix = os.path.join(ICA_AROMAoutDir,'melodic.ica','melodic_mix')
maxRPcorr = BS.feature_time_series(melmix, ICA_inputs.mcImgPar)
melFTmix = os.path.join(ICA_AROMAoutDir,'melodic.ica','melodic_FTmix')
HFC = BS.feature_frequency(melFTmix, float(TR))
motionICs = BS.classification(ICA_AROMAoutDir, maxRPcorr, edgeFract, HFC, csfFract)
BS.denoising(fslDir, ICA_inputs.func, ICA_AROMAoutDir, melmix, denType, motionICs)

denoised_func=os.path.join(ICA_AROMAoutDir,'denoised_func_data_nonaggr.nii.gz')
Temp_Filt_Func = BS.TemporalFilter(denoised_func,Temp_FiltoutDir)

seed='/mnt/cronos/Users/jdkent/VossLabMount/Universal_Software/RestingState2014a/ROIs/LNAcc.nii.gz'
seedname=os.path.basename(seed).replace('.nii.gz','')
Nuis_reg=os.path.join(outDir,'Nuisance_Regression_%s' % (seedname))

 NuisanceReg_func = BS.NuisanceRegression(Temp_Filt_Func,Nrois,ICA_inputs.MNItofunc_warp,Nuis_reg,motion=False)

Prep_BetaSeries=os.path.join(outDir,'Prep_BetaSeries')
if not os.path.isdir(Prep_BetaSeries):
	os.makedirs(Prep_BetaSeries)

EVs=['/home/james/RestingState_dev/FLANKER/custom_timing_files/ev1.txt',
'/home/james/RestingState_dev/FLANKER/custom_timing_files/ev2.txt',
'/home/james/RestingState_dev/FLANKER/custom_timing_files/ev3.txt',
'/home/james/RestingState_dev/FLANKER/custom_timing_files/ev4.txt']

whichEVs=[1, 2, 3]
numrealev=4
tempderiv='-tempderiv'

BS.BetaSeries(Prep_BetaSeries,whichEVs,numrealev,tempderiv=None)
