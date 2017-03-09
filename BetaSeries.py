#!/usr/bin/env python

# Import required modules
import os
import argparse
import commands
import BetaSeries_functions as betafunc
import shutil

# Change to script directory
cwd = os.path.realpath(os.path.curdir)
scriptDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(scriptDir)

#-------------------------------------------- PARSER --------------------------------------------#
#what Do I need as inputs?
#Required
#func
#T1_brain
#T1_head
#EVs (txt file per EV)
#optional
#fieldmaps?

parser = argparse.ArgumentParser(description='Script to run BetaSeries v0.1 beta (\'Beta Series Correlation\') on fMRI data. See the companion manual for further information.')

# Required options                    
reqoptions = parser.add_argument_group('Required arguments')
reqoptions.add_argument('-o', '-out', dest="outDir", required=True, help='Output directory name' )

# Required options in non-Feat mode
nonfeatoptions = parser.add_argument_group('Required arguments - generic mode')
nonfeatoptions.add_argument('-f', '--func',dest="func", required=True, help='Input file name of fMRI data (.nii.gz)')
nonfeatoptions.add_argument('-t', '--t1brain',dest="t1brain",required=True)
nonfeatoptions.add_argument('-T', '--t1head',dest="t1head",required=True)
nonfeatoptions.add_argument('-e', '--EVs',dest="EVs",required=True,nargs='+')
nonfeatoptions.add_argument('-n', '--NROIS',dest="Nrois",required=True)
nonfeatoptions.add_argument('-s', '--Seeds',dest="Seeds",required=True)
nonfeatoptions.add_argument('-w', '--whichEVs',dest="whichEVs",required=True,nargs='+')
nonfeatoptions.add_argument('-tempderiv',dest='tempderiv',action='store_true',
            default=False,help='Include tag if the original design matrix includes temporal derivates. The code assumes that temporal derivatives are immediately after each EV/motion parameter in the Feat design matrix.')

#nrois
#seeds
#don't know if this works
#nonfeatoptions.add_argument('-s', '--smooth',required=False,action='store_false')
#TODO: fieldmap stuff
#nonfeatoptions.add_argument('-mc', dest="mc", required=False, help='motion correction file (e.g., /home/user/PROJECT/SUBJECT.feat/mc/prefiltered_func_data_mcf.par')
#nonfeatoptions.add_argument('-a','-affmat', dest="affmat", default="", help='File name of the mat-file describing the affine registration (e.g., FSL FLIRT) of the functional data to structural space (.mat file). (e.g., /home/user/PROJECT/SUBJECT.feat/reg/example_func2highres.mat')
#nonfeatoptions.add_argument('-w','-warp', dest="warp", default="", help='File name of the warp-file describing the non-linear registration (e.g., FSL FNIRT) of the structural data to MNI152 space (.nii.gz). (e.g., /home/user/PROJECT/SUBJECT.feat/reg/highres2standard_warp.nii.gz')
#nonfeatoptions.add_argument('-m','-mask', dest="mask", default="", help='File name of the mask to be used for MELODIC (denoising will be performed on the original/non-masked input data)')

# Required options in Feat mode
#featoptions = parser.add_argument_group('Required arguments - FEAT mode')
#featoptions.add_argument('-f', '-feat',dest="inFeat", required=False, help='Feat directory name (Feat should have been run without temporal filtering and including registration to MNI152)')

# Optional options
optoptions = parser.add_argument_group('Optional arguments')
optoptions.add_argument('-tr', dest="TR", help='TR in seconds',type=float)
optoptions.add_argument('-den', dest="denType", default="nonaggr", help='Type of denoising strategy: \'no\': only classification, no denoising; \'nonaggr\': non-aggresssive denoising (default); \'aggr\': aggressive denoising; \'both\': both aggressive and non-aggressive denoising (seperately)')
optoptions.add_argument('-md','-meldir', dest="melDir", default="",help='MELODIC directory name, in case MELODIC has been run previously.')
optoptions.add_argument('-dim', dest="dim", default=0,help='Dimensionality reduction into #num dimensions when running MELODIC (default: automatic estimation; i.e. -dim 0)',type=int)
optoptions.add_argument('-eig',dest='eig',action='store_true',default=False)
print '\n------------------------------- RUNNING ICA-AROMA ------------------------------- '
print '--------------- \'ICA-based Automatic Removal Of Motion Artifacts\' --------------- \n'


#--------------------------------------- PARSE ARGUMENTS ---------------------------------------#
args = parser.parse_args()
# Define Variables based on inputs
func=args.func
t1brain=args.t1brain
t1head=args.t1head
EVs=args.EVs
outDir=args.outDir
Nrois_list=args.Nrois
whichEVs=args.whichEVs
eig=args.eig

with open(Nrois_list,'r') as txt:
		Nrois=txt.read().strip().split('\n')

Seeds_list=args.Seeds

with open(Seeds_list,'r') as txt:
		Seeds=txt.read().strip().split('\n')
smooth=True
#make outDir if it doesn't exist
if not os.path.exists(outDir):
		os.makedirs(outDir)
#--------------------------------------- PRE ICA ---------------------------------------#
PreprocessDir=os.path.join(outDir,'preprocess')
ICA_inputs=betafunc.PreICA(PreprocessDir,func,t1brain,t1head,smooth=smooth)



#---------------------------------------- Run ICA-AROMA ----------------------------------------#
##set up variables

fslDir=os.path.join(os.environ["FSLDIR"],'bin','')
ICA_AROMAoutDir=os.path.join(PreprocessDir,"ICA_AROMA")
melDir=os.path.join(ICA_AROMAoutDir,"melodic.ica")

#hacking because of time crunch, manually setting dim and TR
TR=2.0
dim=0
denType='nonaggr'
melIC = os.path.join(ICA_AROMAoutDir,'melodic_IC_thr.nii.gz')
melIC_MNI =  os.path.join(ICA_AROMAoutDir,'melodic_IC_thr_MNI2mm.nii.gz')
#end hack

denoised_func=os.path.join(ICA_AROMAoutDir,'denoised_func_data_nonaggr.nii.gz')
denoised_func_smooth=denoised_func.replace('.nii.gz','_smooth.nii.gz')
denoised_func_nosmooth=denoised_func.replace('.nii.gz','_nosmooth.nii.gz')
if not os.path.isfile(denoised_func_smooth) and not os.path.isfile(denoised_func_nosmooth):
	print 'Step 1) MELODIC'
	if not os.path.isdir(ICA_AROMAoutDir):
		os.makedirs(ICA_AROMAoutDir)
	betafunc.runICA(fslDir, ICA_inputs.func_smooth, ICA_AROMAoutDir,"", ICA_inputs.func_mask, dim, TR)

	print 'Step 2) Automatic classification of the components'
	print '  - registering the spatial maps to MNI'
	betafunc.register2MNI(fslDir, melIC, melIC_MNI, ICA_inputs.functoT1_transform, ICA_inputs.T1toMNI_transform)

	print '  - extracting the CSF & Edge fraction features'
	edgeFract, csfFract = betafunc.feature_spatial(fslDir, ICA_AROMAoutDir, scriptDir, melIC_MNI)

	print '  - extracting the Maximum RP correlation feature'
	melmix = os.path.join(ICA_AROMAoutDir,'melodic.ica','melodic_mix')
	maxRPcorr = betafunc.feature_time_series(melmix, ICA_inputs.mcImgPar)

	print '  - extracting the High-frequency content feature'
	melFTmix = os.path.join(ICA_AROMAoutDir,'melodic.ica','melodic_FTmix')
	HFC = betafunc.feature_frequency(melFTmix, TR)

	print '  - classification'
	motionICs = betafunc.classification(ICA_AROMAoutDir, maxRPcorr, edgeFract, HFC, csfFract)

	if (denType != 'no'):
		print 'Step 3) Data denoising'
		print 'denoising smoothed data'
		betafunc.denoising(fslDir, ICA_inputs.func_smooth, ICA_AROMAoutDir, melmix, denType, motionICs)
		denoised_func=os.path.join(ICA_AROMAoutDir,'denoised_func_data_nonaggr.nii.gz')
		denoised_func_smooth=denoised_func.replace('.nii.gz','_smooth.nii.gz')
		os.rename(denoised_func,denoised_func_smooth)

		print 'denoising unsmoothed data'
		betafunc.denoising(fslDir, ICA_inputs.func_nosmooth, ICA_AROMAoutDir, melmix, denType, motionICs)
		denoised_func_nosmooth=denoised_func.replace('.nii.gz','_nosmooth.nii.gz')
		os.rename(denoised_func,denoised_func_nosmooth)
	# Remove thresholded melodic_IC file
	os.remove(melIC)

	# Revert to old directory
	os.chdir(cwd)

	print '\n----------------------------------- Finished -----------------------------------\n'
#---------------------------------------- Post ICA-AROMA ----------------------------------------#
#denoised_func will be used from here for post ICA processing
#temporally filter the data
Temp_FiltoutDir=os.path.join(PreprocessDir,'Temperal_filter')
if not os.path.isdir(Temp_FiltoutDir):
	os.makedirs(Temp_FiltoutDir)
Temp_Filt_Func = betafunc.TemporalFilter(denoised_func_nosmooth,Temp_FiltoutDir)

#Complete nuisance regression on the data
#Seeds is a list of seed regions
#Nrois multiple txt files with each containing nrois for that seed.
BetaSeriesDir=os.path.join(outDir,'BetaSeries')
for seed,Nroi in zip(Seeds,Nrois):
	seedname=os.path.basename(seed).replace('.nii.gz','')
	if not eig:
		Nuis_reg=os.path.join(PreprocessDir,'Nuisance_Regression_%s' % (seedname))
	elif eig:
		Nuis_reg=os.path.join(PreprocessDir,'Nuisance_Regression_%s_eig' % (seedname))
	else:
		print 'ERROR'
		return 1
	if not os.path.isdir(Nuis_reg):
		os.makedirs(Nuis_reg)
	NuisanceReg_func = betafunc.NuisanceRegression(Temp_Filt_Func,Nroi,ICA_inputs.MNItofunc_warp,Nuis_reg,motion=False,eig)

	#complete betaseries prep on the data
	if not eig:
		Run_BetaSeries=os.path.join(BetaSeriesDir,'BetaSeries_%s' % (seedname))
	elif eig:
		Run_BetaSeries=os.path.join(BetaSeriesDir,'BetaSeries_%s_eig' % (seedname))
	else:
		print 'ERROR'
		
	if not os.path.isdir(Run_BetaSeries):
		os.makedirs(Run_BetaSeries)
	betafunc.MakeModel(NuisanceReg_func,EVs,Run_BetaSeries)

	#hack from time crunch
	whichEVs=[1, 2, 3] #pass in as a option
	numrealev=4 #len(EVs)
	tempderiv='-tempderiv'
	#end hack
	func_data=os.path.join(Run_BetaSeries,'filtered_func_data.nii.gz')
	mask=os.path.join(Run_BetaSeries,'mask.nii.gz')
	if not os.path.islink(func_data):
		os.symlink(NuisanceReg_func,func_data)
	if not os.path.islink(mask):
		os.symlink(ICA_inputs.func_mask,mask)
	#run the betaseries (which expects func and mask in the appropiate place)
	if not os.path.isfile(os.path.join(Run_BetaSeries,'betaseries/ev1_LSS.nii.gz')):
		betafunc.BetaSeries(Run_BetaSeries,whichEVs,numrealev,tempderiv)

	#Extract seed for each ev and run the correlation
	for ev in whichEVs:
		EVLSS=os.path.join(Run_BetaSeries,'betaseries/ev%s_LSS.nii.gz' % (ev))
		Seed_Outdir=os.path.join(outDir,"Results/%s/%s/" % (seedname,"ev"+str(ev)))
		betafunc.SeedCorrelate(EVLSS,seed,Seed_Outdir,ICA_inputs.MNItofunc_warp,ICA_inputs.functoMNI_warp,eig)

#next step. Use seed regions to do correlations

#testcall
#BetaSeries.py -o /home/james/RestingState_dev/FLANKER/outputs/test02282017 -f /home/james/RestingState_dev/FLANKER/inputs/func/pre_sub10_FLANKER_RPI.nii.gz -t /home/james/RestingState_dev/FLANKER/inputs/struct/pre_sub10_MPRAGE-1_RPI_ss.nii.gz -T /home/james/RestingState_dev/FLANKER/inputs/struct/pre_sub10_MPRAGE-1_RPI.nii.gz -e /home/james/RestingState_dev/FLANKER/inputs/custom_timing_files/ev1.txt /home/james/RestingState_dev/FLANKER/inputs/custom_timing_files/ev2.txt /home/james/RestingState_dev/FLANKER/inputs/custom_timing_files/ev3.txt /home/james/RestingState_dev/FLANKER/inputs/custom_timing_files/ev4.txt -n /home/james/RestingState_dev/FLANKER/inputs/nrois/Nrois.txt -s /home/james/RestingState_dev/FLANKER/inputs/seeds/seeds.txt