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
#nrois
#seeds
#don't know if this works
nonfeatoptions.add_argument('-s', '--smooth',required=False,action='store_false')
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
#make outDir if it doesn't exist
if not os.path.exists(outDir):
		os.makedirs(outDir)
#--------------------------------------- PRE ICA ---------------------------------------#
ICA_Inputs=PreICA(outDir,func,t1brain,t1head,smooth=smooth)

# #do motion correction
# mcDir=os.path.join(outDir,'mc')
# if not os.path.exists(mcDir):
# 		os.makedirs(mcDir)

# mcRes=betafunc.motion_correct(func,mcDir):

# #do the registrations
# regDir=os.path.join(outDir,'reg')
# if not os.path.exists(regDir):
# 	os.maskedirs(regDir)

# regRes=betafunc.GenTransforms(regDir,mcRes.mcImgMean,t1brain,t1head)

# #skullstrip the EPI:
# betafunc.skullstripEPI(mcRes.mcImg,t1brain,regRes.T1tofunc_transform)

# #do smoothing if applied (recommended for ICA)
# if smooth:
# 	smoothDir=os.path.join(outDir,'smooth')
# 	betafunc.Spatial_Smooth()


# Define variables based on the type of input (i.e. Feat directory or specific input arguments), and check whether the specified files exist.
# cancel=False
# if args.inFeat:
# 	inFeat = args.inFeat

# 	# Check whether the Feat directory exists
# 	if not os.path.isdir(inFeat): 
# 		print 'The specified Feat directory does not exist.'
# 		print '\n----------------------------- ICA-AROMA IS CANCELED -----------------------------\n'
# 		exit()

# 	# Define the variables which should be located in the Feat directory
# 	inFile = os.path.join(args.inFeat,'filtered_func_data.nii.gz')
# 	mc = os.path.join(args.inFeat,'mc','prefiltered_func_data_mcf.par')
# 	affmat = os.path.join(args.inFeat,'reg','example_func2highres.mat')
# 	warp = os.path.join(args.inFeat,'reg','highres2standard_warp.nii.gz')

# 	# Check whether these files actually exist
# 	if not os.path.isfile(inFile): 
# 		print 'Missing filtered_func_data.nii.gz in Feat directory.'
# 		cancel=True
# 	if not os.path.isfile(mc): 
# 		print 'Missing mc/prefiltered_func_data_mcf.mat in Feat directory.'
# 		cancel=True
# 	if not os.path.isfile(affmat): 
# 		print 'Missing reg/example_func2highres.mat in Feat directory.' 
# 		cancel=True
# 	if not os.path.isfile(warp): 
# 		print 'Missing reg/highres2standard_warp.nii.gz in Feat directory.'
# 		cancel=True
	
# 	# Check whether a melodic.ica directory exists
# 	if os.path.isdir(os.path.join(args.inFeat,'filtered_func_data.ica')):
# 		melDir = os.path.join(args.inFeat,'filtered_func_data.ica')
# 	else: 
# 		melDir=args.melDir
# else:
# 	inFile = args.inFile
# 	mc = args.mc
# 	affmat = args.affmat
# 	warp = args.warp
# 	melDir = args.melDir

# 	# Check whether the files exist
# 	if not inFile:
# 		print 'No input file specified.'
# 	else:
# 		if not os.path.isfile(inFile): 
# 			print 'The specified input file does not exist.'
# 			cancel=True
# 	if not mc:
# 		print 'No mc file specified.'
# 	else:
# 		if not os.path.isfile(mc): 
# 			print 'The specified mc file does does not exist.'
# 			cancel=True
# 	if affmat:
# 		if not os.path.isfile(affmat): 
# 			print 'The specified affmat file does not exist.'
# 			cancel=True
# 	if warp:
# 		if not os.path.isfile(warp): 
# 			print 'The specified warp file does not exist.'
# 			cancel=True

# # Parse the arguments which do not depend on whether a Feat directory has been specified
# outDir = args.outDir
# dim = args.dim
# denType = args.denType

# # Check if the mask exists, when specified.
# if args.mask:
# 	if not os.path.isfile(args.mask):
# 		print 'The specified mask does not exist.'
# 		cancel=True

# # Check if the type of denoising is correctly specified, when specified
# if not (denType == 'nonaggr') and not (denType == 'aggr') and not (denType == 'both') and not (denType == 'no'):
# 	print 'Type of denoising was not correctly specified. Non-aggressive denoising will be run.'
# 	denType='nonaggr'

# # If the criteria for file/directory specifications have not been met. Cancel ICA-AROMA.
# if cancel:
# 	print '\n----------------------------- ICA-AROMA IS CANCELED -----------------------------\n'
# 	exit()

# #------------------------------------------- PREPARE -------------------------------------------#

# # Define the FSL-bin directory
# fslDir = os.path.join(os.environ["FSLDIR"],'bin','')

# # Create output directory if needed
# if not os.path.isdir(outDir):
# 	os.makedirs(outDir)

# # Get TR of the fMRI data, if not specified
# if args.TR:
# 	TR = args.TR
# else:
# 	cmd = ' '.join([os.path.join(fslDir,'fslinfo'), 
# 		inFile, 
# 		'| grep pixdim4 | awk \'{print $2}\''])
# 	TR=float(commands.getoutput(cmd))

# # Check TR
# if TR == 1:
# 	print 'Warning! Please check whether the determined TR (of ' + str(TR) + 's) is correct!\n'
# elif TR == 0:
# 	print 'TR is zero. ICA-AROMA requires a valid TR and will therefore exit. Please check the header, or define the TR as an additional argument.\n----------------------------- ICA-AROMA IS CANCELED -----------------------------\n'

# # Define/create mask. Either by making a copy of the specified mask, or by creating a new one.
# mask = os.path.join(outDir,'mask.nii.gz')
# if args.mask:
# 	shutil.copyfile(args.mask, mask)
# else:
# 	# If a Feat directory is specified, and an example_func is present use example_func to create a mask
# 	if args.inFeat and os.path.isfile(os.path.join(inFeat,'example_func.nii.gz')):
# 		os.system(' '.join([os.path.join(fslDir,'bet'),
# 			os.path.join(inFeat,'example_func.nii.gz'), 
# 			os.path.join(outDir,'bet'),
# 			'-f 0.3 -n -m -R']))
# 		os.system(' '.join(['mv',
# 			os.path.join(outDir,'bet_mask.nii.gz'), 
# 			mask]))
# 		if os.path.isfile(os.path.join(outDir,'bet.nii.gz')):
# 			os.remove(os.path.join(outDir,'bet.nii.gz'))
# 	else:
# 		if args.inFeat:
# 			print ' - No example_func was found in the Feat directory. A mask will be created including all voxels with varying intensity over time in the fMRI data. Please check!\n'
# 		os.system(' '.join([os.path.join(fslDir,'fslmaths'),
# 			inFile,
# 			'-Tstd -bin',
# 			mask]))


#---------------------------------------- Run ICA-AROMA ----------------------------------------#
##set up variables

fslDir=os.path.join(os.environ["FSLDIR"],'bin','')
ICA_AROMAoutDir=os.path.join(outDir,"ICA_AROMA")
melDir=os.path.join(ICA_AROMAoutDir,"melodic.ica")

#hacking because of time crunch, manually setting dim and TR
TR=2.0
dim=0
denType='nonaggr'
melIC = os.path.join(ICA_AROMAoutDir,'melodic_IC_thr.nii.gz')
melIC_MNI =  os.path.join(ICA_AROMAoutDir,'melodic_IC_thr_MNI2mm.nii.gz')
#end hack


if smooth:
	print 'Step 1) MELODIC'
	if not os.path.isdir(ICA_AROMAoutDir):
		os.makedirs(ICA_AROMAoutDir)
	betafunc.runICA(fslDir, ICA_inputs.func, ICA_AROMAoutDir, melDir, ICA_inputs.func_mask, dim, TR)

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
		betafunc.denoising(fslDir, ICA_inputs.func, ICA_AROMAoutDir, melmix, denType, motionICs)

	# Remove thresholded melodic_IC file
	os.remove(melIC)

	# Revert to old directory
	os.chdir(cwd)

	denoised_func=os.path.join(ICA_AROMAoutDir,'denoised_func_data_nonaggr.nii.gz')
	print '\n----------------------------------- Finished -----------------------------------\n'
else:

	#no smoothing was ran, so no ICA_AROMA was ran, presumed it was already ran
	#extract the IC's that are associated with motion
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
	ICA_AROMA_ns_outDir=os.path.join(outDir,"ICA_AROMA_no_smooth")
	print 'Step 3) Data denoising'
	betafunc.denoising(fslDir, ICA_inputs.func, ICA_AROMA_ns_outDir, melmix, denType, motionICs)

	denoised_func=os.path.join(ICA_AROMA_ns_outDir,'denoised_func_data_nonaggr.nii.gz')

	print '\n----------------------------------- Finished -----------------------------------\n'
#---------------------------------------- Post ICA-AROMA ----------------------------------------#
#denoised_func will be used from here for post ICA processing
#temporally filter the data
Temp_FiltoutDir=os.path.join(outDir,'Temperal_filter')
if not os.path.isdir(Temp_FiltoutDir):
	os.makedirs(Temp_FiltoutDir)
Temp_Filt_Func = betafunc.TemporalFilter(denoised_func,Temp_FiltoutDir)

#Complete nuisance regression on the data
seedname=os.path.basename(seed).replace('.nii.gz','')
Nuis_reg=os.path.join(outDir,'Nuisance_Regression_%s' % (seedname))
if not os.path.isdir(Nuis_reg):
	os.makedirs(Nuis_reg)
NuisanceReg_func = betafunc.NuisanceRegression(Temp_Filt_Func,nrois,ICA_Inputs.MNItofunc_warp,Nuis_reg,motion=False)

#complete betaseries prep on the data
Prep_BetaSeries=os.path.join(outDir,'Prep_BetaSeries')
if not os.path.isdir(Prep_BetaSeries):
	os.makedirs(Prep_BetaSeries)
betafunc.MakeModel(NuisanceReg_func,EVs,Prep_BetaSeries)

#run betaseries
#hack from time crunch
whichEVs=[1, 2, 3]
numrealev=4
tempderiv='-tempderiv'
#end hack

betafunc.BetaSeries(Prep_BetaSeries,whichEVs,numrealev,tempderiv=None)


#next step. Use seed regions to do correlations
