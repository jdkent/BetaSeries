#!/usr/bin/env python
def getmiddlevolume(func):
    from nibabel import load
    funcfile = func
    if isinstance(func, list):
        funcfile = func[0]
    _, _, _, timepoints = load(funcfile).shape
    return int(timepoints / 2) - 1


def makeMask(image):
	#binarizes a 3d volume
	import os
	import collections
	import nipype.interfaces.fsl as fsl
	import string
	#Name the output (NIFTI)
	currentdir=os.getcwd()
	outdir=os.path.dirname(image)

	os.chdir(outdir)

	MaskFile=string.replace(image,".nii.gz","_mask.nii.gz")

	#set up the command
	MaskCmd=fsl.ImageMaths()
	MaskCmd.inputs.in_file=image
	MaskCmd.inputs.out_file=MaskFile
	MaskCmd.inputs.op_string='-bin'

	try:
		MaskCmd.run()
	except:
		os.chdir(currentdir)
		print 'ERROR: ' + MaskCmd.cmdline + '\n'
		return 1

	print 'SUCCESS: ' + MaskCmd.cmdline + '\n'
	os.chdir(currentdir)
	return MaskFile


def MotionCorrect(func,outDir):
	import os
	import subprocess
	from collections import namedtuple
	try:
		import nipype.interfaces.fsl as fsl
		import nipype.interfaces.afni as afni
		import numpy as np
	except:
		print "Have you installed nipype? 'sudo pip install nipype'"
		return 1

	mcImg=os.path.join(outDir,'mcImg.nii.gz')
	mcImgRaw=os.path.join(outDir,'mcImg_raw.par')
	mcImgMat=os.path.join(outDir,'mcImg_mat.aff12.1D')
	mcImg1D=os.path.join(outDir,'mcImg1D.mat')
	mcImgMaxDisp=os.path.join(outDir,'mcImgMaxDisp.mat')
	mcImgMean=os.path.join(outDir,'mcImgmean.nii.gz')
	mcImgDegPar=os.path.join(outDir,'mcImg_deg.par')
	mcImgPar=os.path.join(outDir,'mcImg.par')
	mcImgMMPar=os.path.join(outDir,'mcImg_mm.par')
	mcImgABS=os.path.join(outDir,'mcImg_abs.rms')
	mcImgDerivPar=os.path.join(outDir,'mcImg_deriv.par')
	mcImgRel=os.path.join(outDir,'mcImg_rel.rms')
	rotPNG=os.path.join(outDir,'rot.png')
	transPNG=os.path.join(outDir,'trans.png')
	rotMMPNG=os.path.join(outDir,'rot_mm.png')
	rotTRANSPNG=os.path.join(outDir,'rot_trans.png')
	ABSPNG=os.path.join(outDir,'abs_mot.png')
	RELPNG=os.path.join(outDir,'rel_mot.png')
	DISPNG=os.path.join(outDir,'disp.png')
	############# Set up OutDir and Log file ###########
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	currentdir=os.getcwd()
	os.chdir(outDir)
	#################################################################
	#logfile
	MCLogFile=os.path.join(outDir,'mc.log')

	if os.path.isfile(MCLogFile):
		MCLog = open(MCLogFile,"a")
		MCLog.write("-------------New Process-------------\n")
	else:
		MCLog = open(MCLogFile,"w")

	if os.path.isfile(mcImg):
		msg="- mcImg.nii.gz exists, skipping this step\n"
		print msg
		MCLog.write(msg)
	else:
		#get the middle volume
		try:
			halfpoint=getmiddlevolume(func)
		except:
			msg='ERROR: getmiddlevolume\n'
			print msg
			MCLog.write(msg)
		#complete motion correction.
		try:
			MCcmd=afni.Volreg()
			MCcmd.inputs.in_file=func
			MCcmd.inputs.oned_file=mcImg1D
			MCcmd.inputs.md1d_file=mcImgMaxDisp
			MCcmd.inputs.oned_matrix_save=mcImgMat
			MCcmd.inputs.out_file=mcImg
			MCcmd.inputs.outputtype='NIFTI_GZ'
			MCcmd.inputs.args='-Fourier -twopass -base %d -dfile %s' % (halfpoint, mcImgRaw)
			MCcmd.inputs.timeshift=False
			MCcmd.inputs.zpad=4
			
			#run command
			MCcmd_Result=MCcmd.run()
		except:
			print 'ERROR: ' + MCcmd.cmdline + "\n"
			MCLog.write('ERROR: ' + MCcmd.cmdline + "\n")
			MCLog.close()
			return 1

		#log that flirt finished successfully			
		MCLog.write('SUCCESS: ' + MCcmd.cmdline + "\n")

	#Take the mean of mcImg.nii.gz (for use in reg)
	if os.path.isfile(mcImgMean):
		msg="- mcImgMean.nii.gz exists, skipping this step\n"
		print msg
		MCLog.write(msg)
	else:
		 #get the mean from mcImg.nii.gz
		 try:
		 	 MeanCmd=fsl.MeanImage()
		 	 MeanCmd.inputs.in_file=mcImg
		 	 MeanCmd.inputs.dimension='T'
		 	 MeanCmd.inputs.out_file=mcImgMean
		 	 MeanCmd.inputs.nan2zeros=True
		 	 MeanCmd.inputs.internal_datatype='float'
		 	 MeanCmd.inputs.output_datatype='float'
		 	 MeanCmd.run()

		 except:
		 	 msg='ERROR: ' + MeanCmd.cmdline +'\n'
		 	 print msg
		 	 MCLog.write(msg)
		 	 MCLog.close()
		 	 return 1

		 #log that flirt finished successfully			
		 MCLog.write('SUCCESS: ' + MeanCmd.cmdline + "\n")



	#read in mcImg_raw.par
	#Save out mcImg.par (like fsl) with only the translations and rotations
    #mcflirt appears to have a different rotation/translation order.  Reorder 3dvolreg output to match "RPI" FSL ordering
    ##AFNI ordering
    #roll  = rotation about the I-S axis }
 	#pitch = rotation about the R-L axis } degrees CCW
  	#yaw   = rotation about the A-P axis }
  	#dS  = displacement in the Superior direction  }
  	#dL  = displacement in the Left direction      } mm
  	#dP  = displacement in the Posterior direction }
	mcImgRawFile=open(mcImgRaw,'r')
	mcImgRaw_np=np.array([[float(x) for x in line.strip().split()] for line in mcImgRawFile])
	mcImgRawFile.close()

	mcImgDeg_np=mcImgRaw_np[:,[2,3,1,5,6,4]]
	mcImgDegPar_id=open(mcImgDegPar, 'w')
	np.savetxt(mcImgDegPar_id,mcImgDeg_np,fmt='%.10f')
	mcImgDegPar_id.close()

	#Need to convert rotational parameters from degrees to radians
    #rotRad= (rotDeg*pi)/180
    #pi=3.14159
	pi=3.14159
	trans_mcImgPar_np=np.multiply(mcImgDeg_np[:,[0,1,2]],(pi/180))
	notrans_mcImgPar_np=mcImgDeg_np[:,[3,4,5]]
	mcImgPar_np=np.concatenate((trans_mcImgPar_np,notrans_mcImgPar_np),axis=1)
	mcImgPar_id=open(mcImgPar, 'w')
	np.savetxt(mcImgPar_id,mcImgPar_np,fmt='%.10f')
	mcImgPar_id.close()

	#Need to create a version where ALL (rotations and translations) measurements are in mm.  Going by Power 2012 Neuroimage paper, radius of 50mm.
    #Convert degrees to mm, leave translations alone.
    #rotDeg= ((2r*Pi)/360) * Degrees = Distance (mm)
    #d=2r=2*50=100
    #pi=3.14159
	d=100
	trans_mcImgMM_np=np.multiply(mcImgDeg_np[:,[0,1,2]],(d*pi/360))
	notrans_mcImgMM_np=mcImgDeg_np[:,[3,4,5]]
	mcImgMM_np=np.concatenate((trans_mcImgMM_np,notrans_mcImgMM_np),axis=1)
	mcImgMM_id=open(mcImgMMPar, 'w')
	np.savetxt(mcImgMM_id,mcImgMM_np,fmt='%.10f')
	mcImgMM_id.close()

	#Cut motion parameter file into 6 distinct TR parameter files
	for x in range(0,6):
		index=x+1
		mc_file=(os.path.join(outDir,"mc" + str(index) + ".txt"))
		mc_file_id=open(mc_file, 'w')
		np.savetxt(mc_file_id,mcImgPar_np[:,x],fmt='%.10f')
		mc_file_id.close()


	#Absolute Displacement
	subprocess.check_output("awk '{print (sqrt(0.2*80^2*((cos($1)-1)^2+(sin($1))^2 + (cos($2)-1)^2 + (sin($2))^2 + (cos($3)-1)^2 + (sin($3)^2)) + $4^2+$5^2+$6^2))}' %s > %s" % (mcImgPar, mcImgABS),stderr=subprocess.STDOUT, shell=True)

	
    #Create the relative displacement .par file from the input using AFNI's 1d_tool.py to first calculate the derivatives
	subprocess.check_output("1d_tool.py -infile %s -set_nruns 1 -derivative -write %s -overwrite" % (mcImgPar, mcImgDerivPar),stderr=subprocess.STDOUT, shell=True)

	#Relative Displacement
	subprocess.check_output("awk '{print (sqrt(0.2*80^2*((cos($1)-1)^2+(sin($1))^2 + (cos($2)-1)^2 + (sin($2))^2 + (cos($3)-1)^2 + (sin($3)^2)) + $4^2+$5^2+$6^2))}' %s > %s" % (mcImgDerivPar, mcImgRel),stderr=subprocess.STDOUT, shell=True)

	#rotation plot (degrees)
	RotPlot=fsl.PlotMotionParams()
	RotPlot.inputs.in_file=mcImgPar
	RotPlot.inputs.in_source='fsl'
	RotPlot.inputs.plot_type='rotations'
	RotPlot.inputs.plot_size=(300,800)
	RotPlot.inputs.out_file=rotPNG
	RotPlot_Result=RotPlot.run()

	#translation plot (mm)
	TransPlot=fsl.PlotMotionParams()
	TransPlot.inputs.in_file=mcImgPar
	TransPlot.inputs.in_source='fsl'
	TransPlot.inputs.plot_type='translations'
	TransPlot.inputs.plot_size=(300,800)
	TransPlot.inputs.out_file=transPNG
	TransPlot_Result=TransPlot.run()

	#rotation plot (mm)
	RotMMPlot=fsl.PlotMotionParams()
	RotMMPlot.inputs.in_file=mcImgMMPar
	RotMMPlot.inputs.in_source='fsl'
	RotMMPlot.inputs.plot_type='rotations'
	RotMMPlot.inputs.plot_size=(300,800)
	RotMMPlot.inputs.out_file=rotMMPNG
	RotMMPlot_Result=RotMMPlot.run()

	#displacement (mm)
	RotTransMMPlot=fsl.PlotMotionParams()
	RotTransMMPlot.inputs.in_file=mcImgMMPar
	RotTransMMPlot.inputs.in_source='fsl'
	RotTransMMPlot.inputs.plot_type='displacement'
	RotTransMMPlot.inputs.plot_size=(300,800)
	RotTransMMPlot.inputs.out_file=rotMMPNG
	RotTransMMPlot_Result=RotMMPlot.run()
	
	#absolute motion
	AbsPlot=fsl.PlotMotionParams()
	#AbsRelPlot.inputs.in_file=mcImgABS
	AbsPlot.inputs.in_file=mcImgABS
	AbsPlot.inputs.in_source='fsl'
	#AbsRelPlot.inputs.args='-i %s,%s -a absolute,relative -t "3dvolreg estimated mean displacement (mm)"' % (mcImgABS, mcImgRel)
	AbsPlot.inputs.plot_type='displacement'
	AbsPlot.inputs.plot_size=(300,800)
	AbsPlot.inputs.out_file=ABSPNG
	AbsPlot_Result=AbsPlot.run()

	#relative motion
	RelPlot=fsl.PlotMotionParams()
	RelPlot.inputs.in_file=mcImgRel
	RelPlot.inputs.in_source='fsl'
	RelPlot.inputs.plot_type='displacement'
	RelPlot.inputs.plot_size=(300,800)
	RelPlot.inputs.out_file=RELPNG
	RelPlot_Result=RelPlot.run()

	#write finished imformation to log and close log.
	msg="SUCCESS: MotionCorrect"
	print msg
	MCLog.write(msg)
	MCLog.close()

	os.chdir(currentdir)

	#return success
	outputs = namedtuple("outfiles", ["mcImg","mcImgMean","mcImgPar"])
	return outputs(mcImg,mcImgMean,mcImgPar)

# Functions for ICA-AROMA v0.3 beta
def GenTransforms(outDir,ex_func,T1_brain,T1_head,fmap=None,fmapmag=None,fmapmagbrain=None,echospacing=None,pedir=None):
	"""################
	   #function: GenTransforms
	   ################
	   #purpose: creates the transforms to move func to MNI (and vice versa)
	   ################
	   #preconditions: fsl is installed
	   ################
	   #input: fslDir- full path to the overview of FSL (e.g. /usr/local/fsl)
	   #       outDir- full path to the directory where the transforms are made 
	   #       ex_func- full path to an example 3D functional image
	   #       T1_brain- full path to the T1_brain
	   #       T1_head- full path to the T1_head
	   # Optional arguments (for fieldmap processing)
	   #	   fmap - full path to fmap
	   #	   fmapmag - full path to field map magnitude (head)
	   #	   fmapmagbrain - full path to field map magnitude (brain)
	   #	   echospacing - dwell time
	   #	   pedir - phase encoding direction
	   ################
	   #output: T1toMNI_warp.nii.gz
	   #        T1toMNI.mat
	   #		T1toMNI_flirt
	   #		T1toMNI_fnirt
	   #		T1toMNIfnirt.log
	   #        functoT1.mat
	   #        functoMNI_warp.nii.gz
	   #        MNItoT1_warp.nii.gz
	   #        MNItoT1.mat
	   #        MNItofunc_warp.nii.gz
	   #		T1tofunc.mat
	   #		OPTIONAL (fieldmap)
	   #		T1tofunc_warp
	   #		functoT1_warp
	   ################"""

	# Import needed modules
	import os
	import string
	from collections import namedtuple
	try:
		import nipype.interfaces.fsl as fsl
	except:
		print "Have you installed nipype? 'sudo pip install nipype'"
		return 1

	############# Set up OutDir and Log file ###########
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	currentdir=os.getcwd()
	os.chdir(outDir)
	#################################################################
	#logfile
	RegLogFile=os.path.join(outDir,'reg.log')

	if os.path.isfile(RegLogFile):
		RegLog = open(RegLogFile,"a")
		RegLog.write("-------------New Process-------------\n")
	else:
		RegLog = open(RegLogFile,"w")
	################# Output Variable Names ###################
	#define fsl pathnames
	#fsl bin (ubuntu installation)
	#fslbin=os.path.join(fslDir,'bin/')

	#fsl MNI files (ubuntu installation)
	#fslstandards=os.path.join(fslDir,'data/standard/')
	fslstandards=os.path.join(os.environ['FSLDIR'],'data/standard')

	# Define output files
	#T1 to MNI
	T1toMNI_warp=os.path.join(outDir,'T1toMNI_warp.nii.gz')
	T1toMNI_mat=os.path.join(outDir,'T1toMNI.mat')
	T1toMNI_flirt=os.path.join(outDir,'T1toMNI_flirt.nii.gz')
	T1toMNI_fnirt=os.path.join(outDir,'T1toMNI_fnirt.nii.gz')
	T1toMNI_fnirt_log=os.path.join(outDir,'T1toMNIfnirt.log')
	#MNI to T1
	MNItoT1_warp=os.path.join(outDir,'MNItoT1_warp.nii.gz')
	MNItoT1_mat=os.path.join(outDir,'MNItoT1.mat')

	#func to T1
	functoT1_base=os.path.join(outDir,'functoT1')
	functoT1_mat=os.path.join(outDir,'functoT1.mat')

	#T1 to func
	T1tofunc_mat=os.path.join(outDir,'T1tofunc.mat')

	#func to MNI
	functoMNI_warp=os.path.join(outDir,'functoMNI_warp.nii.gz')

	#MNI to func
	MNItofunc_warp=os.path.join(outDir,'MNItofunc_warp.nii.gz')

	
	
	#### Check if there is special fieldmap processing ####
	if not fmap==None and not fmapmag==None and not fmapmagbrain==None and not echospacing==None and not pedir==None:
		T1tofunc_warp=os.path.join(outDir,'T1tofunc_warp.nii.gz')
		functoT1_warp=os.path.join(outDir,'functoT1_warp.nii.gz')
		fieldmap=True
	else:
		fieldmap=False

	RegLog.write("Fieldmap: " + str(fieldmap) + "\n")


	#T1toMNI transforms
	if os.path.isfile(T1toMNI_warp) and os.path.isfile(T1toMNI_mat):
			 msg="- Warp and mat from T1 to MNI already exists, skipping this step\n"
			 print msg
			 RegLog.write(msg)
	else: 
		#flirt (mat) linear transform
		try:
			T1toMNIflt = fsl.FLIRT()
			T1toMNIflt.inputs.in_file=T1_brain
			T1toMNIflt.inputs.reference=os.path.join(fslstandards,'MNI152_T1_2mm_brain.nii.gz')
			T1toMNIflt.inputs.out_matrix_file=T1toMNI_mat
			T1toMNIflt.inputs.out_file=T1toMNI_flirt
			
			#print the command being ran:
			#print T1toMNIflt.cmdline

			#run the command:
			T1toMNIflt_Result = T1toMNIflt.run()
		except:
			RegLog.write('ERROR: ' + T1toMNIflt.cmdline + "\n")
			RegLog.close()
			return 1
		#log that flirt finished successfully			
		RegLog.write('SUCCESS: ' + T1toMNIflt.cmdline + "\n")


		#fnirt (warp) nonlinear transform
		try:
			T1toMNIfnt=fsl.FNIRT()
			T1toMNIfnt.inputs.in_file=T1_head
			T1toMNIfnt.inputs.ref_file=os.path.join(fslstandards,'MNI152_T1_2mm.nii.gz')
			T1toMNIfnt.inputs.affine_file=T1toMNI_mat
			T1toMNIfnt.inputs.config_file='T1_2_MNI152_2mm'
			T1toMNIfnt.inputs.fieldcoeff_file=T1toMNI_warp
			T1toMNIfnt.inputs.warped_file=T1toMNI_fnirt
			T1toMNIfnt.inputs.jacobian_range=(.1,10)
			T1toMNIfnt.inputs.log_file=T1toMNI_fnirt_log
			#print the command:
			#print T1toMNIfnt.cmdline

			#run the command
			T1toMNIfnt_Result = T1toMNIfnt.run()
		except:
			RegLog.write('ERROR: ' + T1toMNIfnt.cmdline + "\n")
			RegLog.close()
			return 1
		#log that fnirt finished successfully
		RegLog.write('SUCCESS: ' + T1toMNIfnt.cmdline + "\n")

	#MNItoT1 transforms
	if os.path.isfile(MNItoT1_warp) and os.path.isfile(MNItoT1_mat):
		msg="- Warp and mat from MNI to T1 already exists, skipping this step\n"
		print msg
		RegLog.write(msg)
	else:
		#inverse the T1toMNI linear transform
		try:
			MNItoT1xfm=fsl.ConvertXFM()
			MNItoT1xfm.inputs.in_file=T1toMNI_mat
			MNItoT1xfm.inputs.invert_xfm=True
			MNItoT1xfm.inputs.out_file=MNItoT1_mat

			#print the command:
			#print MNItoT1xfm.cmdline

			#run the command
			MNItoT1xfm_Result = MNItoT1xfm.run()
		except:
			RegLog.write('ERROR: ' + MNItoT1xfm.cmdline + "\n")
			RegLog.close()
			return 1
		#log that fnirt finished successfully
		RegLog.write('SUCCESS: ' + MNItoT1xfm.cmdline + "\n")
		
		try:
			MNItoT1invw=fsl.InvWarp()
			MNItoT1invw.inputs.warp=T1toMNI_warp
			MNItoT1invw.inputs.reference=T1_head
			MNItoT1invw.inputs.inverse_warp=MNItoT1_warp

			#print the command:
			#print MNItoT1invw.cmdline

			#run the command
			MNItoT1invw_Result = MNItoT1invw.run()
		except:
			RegLog.write('ERROR: ' + MNItoT1invw.cmdline + "\n")
			RegLog.close()
			return 1
		#log that fnirt finished successfully
		RegLog.write('SUCCESS: ' + MNItoT1invw.cmdline + "\n")
	

	#functoT1 transform
	if os.path.isfile(functoT1_mat):
		if fieldmap and os.path.isfile(functoT1_warp):
			msg='- Warp and mat from func to T1 already exists, skipping this step\n'
		else:
			msg='- Mat from func to T1 already exists, skipping this step\n'
		print msg
		RegLog.write(msg)
	else:
		#start making the epireg command
		functoT1epireg=fsl.EpiReg()

		#special fieldmap options
		if fieldmap:
			functoT1epireg.inputs.fmap=fmap
			functoT1epireg.inputs.fmapmag=fmapmag
			functoT1epireg.inputs.fmapmagbrain=fmapmagbrain
			functoT1epireg.inputs.echospacing=echospacing
			functoT1epireg.inputs.pedir=pedir
		#default arguments
		functoT1epireg.inputs.epi=ex_func
		functoT1epireg.inputs.t1_head=T1_head
		functoT1epireg.inputs.t1_brain=T1_brain
		functoT1epireg.inputs.out_base=functoT1_base

		try:
			functoT1epireg.run()
		except:
			RegLog.write('ERROR: ' + functoT1epireg.cmdline + "\n")
			RegLog.close()
			return 1
		#log that epi_reg finished successfully
		RegLog.write('SUCCESS: ' + functoT1epireg.cmdline + "\n")

	#T1tofunc transform
	if os.path.isfile(T1tofunc_mat):
		if fieldmap and os.path.isfile(T1tofunc_warp):
			msg='- Warp and mat from T1 to func already exists, skipping this step\n'
		else:
			msg='- Mat from T1 to func already exists, skipping this step\n'
		print msg
		RegLog.write(msg)

	else:
		#do invwarp if there is a fieldmap
		if fieldmap:
			T1tofuncinvw=fsl.InvWarp()
			T1tofuncinvw.inputs.warp=functoT1_warp
			T1tofuncinvw.inputs.reference=ex_func
			T1tofuncinvw.inputs.inverse_warp=T1tofunc_warp

			try:
				T1tofuncinvw.run()
			except:
				RegLog.write('ERROR: ' + T1tofuncinvw.cmdline + "\n")
				RegLog.close()
				return 1

			#log that invwarp finished successfully
			RegLog.write('SUCCESS: ' + T1tofuncinvw.cmdline + "\n")

		#convert xfm for the mat file
		T1tofuncxfm=fsl.ConvertXFM()
		T1tofuncxfm.inputs.in_file=functoT1_mat
		T1tofuncxfm.inputs.invert_xfm=True
		T1tofuncxfm.inputs.out_file=T1tofunc_mat

		try:
			T1tofuncxfm.run()
		except:
			RegLog.write('ERROR: ' + T1tofuncxfm.cmdline + "\n")
			RegLog.close()
			return 1

		#log that invwarp finished successfully
		RegLog.write('SUCCESS: ' +  T1tofuncxfm.cmdline + "\n")


	#sum warps/transforms to get MNItoEPI warp
	if os.path.isfile(MNItofunc_warp):
		msg='- Warp from MNI to EPI already exists, skipping this step\n'
		print msg
		RegLog.write(msg)

	else:

		#setup command
		MNItofuncWarp=fsl.ConvertWarp()
		#special fieldmap option
		if fieldmap:
			MNItofuncWarp.inputs.warp2=T1tofunc_warp
		else:
			MNItofuncWarp.inputs.postmat=T1tofunc_mat

		#default
		MNItofuncWarp.inputs.reference=ex_func
		MNItofuncWarp.inputs.warp1=MNItoT1_warp
		MNItofuncWarp.inputs.out_relwarp=True
		MNItofuncWarp.inputs.out_file=MNItofunc_warp

		try:
			MNItofuncWarp.run()
		except:
			RegLog.write('ERROR: ' + MNItofuncWarp.cmdline + "\n")
			RegLog.close()
			return 1

		#log that ConvertWarp finished successfully
		RegLog.write('SUCCESS: ' +  MNItofuncWarp.cmdline + "\n")

	#invert the MNItofunc warp to make the functoMNI warp
	if os.path.isfile(functoMNI_warp):
		msg="- Warp from func to MNI already exists, skipping this step\n"
		print msg
		RegLog.write(msg)
	else:
		functoMNIinvw=fsl.InvWarp()
		functoMNIinvw.inputs.warp=MNItofunc_warp
		functoMNIinvw.inputs.reference=os.path.join(fslstandards,'MNI152_T1_2mm.nii.gz')
		functoMNIinvw.inputs.inverse_warp=functoMNI_warp

		try:
			functoMNIinvw.run()
		except:
			RegLog.write('ERROR: ' + functoMNIinvw.cmdline + "\n")
			RegLog.close()
			return 1

		#log that invwarp finished successfully
		RegLog.write('SUCCESS: ' + functoMNIinvw.cmdline + "\n")

	#finish the function	
	msg='SUCCESS: GetTransforms\n'
	os.chdir(currentdir)
	print msg
	RegLog.write(msg)

	os.chdir(currentdir)

	outputs = namedtuple("outfiles", ["T1tofunc_transform","functoT1_transform","T1toMNI_transform","MNItoT1_transform","MNItofunc_warp","functoMNI_warp"])
	if fieldmap:
		return outputs(T1tofunc_warp,functoT1_warp,T1toMNI_warp,MNItoT1_warp,MNItofunc_warp,functoMNI_warp)
	else:
		return outputs(T1tofunc_mat,functoT1_mat,T1toMNI_warp,MNItoT1_warp,MNItofunc_warp,functoMNI_warp)

#idk about the inputs, I think this way minimizes the processing?
def Spatial_Smooth(func_brain,func_mask,func_mean,outDir):
	import nipype.interfaces.fsl as fsl
	import string
	from collections import namedtuple
	import os
	#currentdir
	currentdir=os.getcwd()
	#outdir
	os.chdir(outDir)
	#func is motion corrected and stripped.	
	#mean=$(fslstats ${mean_func} -k ${mask} -p 50)
	#brightness threshold=mean*.75
	#susan $1 ${brightness_threshold} 2.5 3 1 1 ${mean_func} ${brightness_threshold} ${smoothName}
	#http://nipype.readthedocs.io/en/latest/interfaces/generated/nipype.interfaces.fsl.preprocess.html
	smooth_out=os.path.join(outDir,"func_smooth.nii.gz")
	meanCmd=fsl.ImageStats(in_file=func_mean,op_string='-k %s -p 50' % func_mask)

	#derive brightness threshold (determined by featlib.tcl)
	try:	
		threshold=meanCmd.run().outputs.out_stat*.75
	except:
		print '- ERROR, could not determine threshold for spatial smoothing'
		return 1

	#run spatial smoothing on the data using fsl's SUSAN
   	spatialSmoothCmd=fsl.SUSAN(brightness_threshold=threshold,fwhm=2.5,in_file=func_brain,dimension=3,usans=[(func_mean,threshold)],out_file=smooth_out)
   	try:
   		spatialSmoothCmd.run()
   	except:
   		print '- ERROR, could not run ' + spatialSmoothCmd.cmdline

   	#print success
   	print '- SUCCESS Spatial_Smooth: %s\n%s' % (meanCmd.cmdline, spatialSmoothCmd.cmdline)
 	os.chdir(currentdir)
   	outputs = namedtuple("outfiles", ["smoothed_output"])
   	return outputs(smooth_out)

#test: BetaSeries_functions.skullstripEPI('/home/james/RestingState_dev/test/mc/mcImg.nii.gz','/home/james/RestingState_dev/data/pre_sub10_MPRAGE-1_RPI_ss.nii.gz','/home/james/RestingState_dev/reg/T1tofunc.mat')
def skullstripEPI(func,T1_brain=None,T1tofunc_transform=None,EPI_mask=None):
	#T1_mask (Assumes NIFTI)
	import os
	import string
	from collections import namedtuple
	import nipype.interfaces.fsl as fsl

	#current directory
	currentdir=os.getcwd()
	#outdir
	outDir=os.path.dirname(func)
	#change directory
	os.chdir(outDir)

	#check transform
	if T1tofunc_transform is not None:
		transform_type=os.path.basename(T1tofunc_transform)
	else:
		transform_type=None
		
	if transform_type is not None and '.mat' in transform_type:
		T1tofunc_mat=T1tofunc_transform
		T1tofunc_warp=None
	elif transform_type is not None and '.nii.gz' in transform_type:
		T1tofunc_mat=None
		T1tofunc_warp=T1tofunc_transform
	else:
		T1tofunc_mat=None
		T1tofunc_warp=None
	#T1_mask=string.replace(T1_brain,".nii.gz","_mask.nii.gz")
	#in case the user doesn'tsupplies a mask
	if EPI_mask is None:
		EPI_mask=string.replace(func,".nii.gz","_mask.nii.gz")
		#make sure necessary inputs are filled in
		if T1_brain is None or (T1tofunc_mat is None and T1tofunc_warp is None):
			print '- ERROR, either specify EPI_mask or T1_brain and (T1tofunc_mat or T1tofunc_warp)'
			return 1
	#output image
	func_brain=string.replace(func,".nii.gz","_brain.nii.gz")
	if T1_brain is not None:
		T1_mask=makeMask(T1_brain)

	if T1tofunc_warp is not None and os.path.isfile(T1tofunc_warp):
		T1tofuncMask=fsl.ApplyWarp()
		T1tofuncMask.inputs.in_file=T1_mask
		T1tofuncMask.inputs.out_file=EPI_mask
		T1tofuncMask.inputs.ref_file=func
		T1tofuncMask.inputs.field_file=T1tofunc_warp
		T1tofuncMask.inputs.interp='nn'
		T1tofuncMask.inputs.datatype='char'

		try:
			T1tofuncMask.run()
		except:
			print 'ERROR: ' + T1tofuncMask.cmdline + "\n"
			return 1
		#log that ApplyWarp finished successfully
		print 'SUCCESS: ' +  T1tofuncMask.cmdline + "\n"
	else:
  		 
  		#tApplyXfm will be deprecated, eventually make it fsl.ApplyXFM
  		if not os.path.isfile(EPI_mask):
	  		T1tofuncMask=fsl.ApplyXfm()
	  		T1tofuncMask.inputs.in_file=T1_mask
			T1tofuncMask.inputs.out_file=EPI_mask
			T1tofuncMask.inputs.reference=func
			T1tofuncMask.inputs.in_matrix_file=T1tofunc_mat
			T1tofuncMask.inputs.interp='nearestneighbour'
			T1tofuncMask.inputs.datatype='char'

			try:
				T1tofuncMask.run()
			except:
				print 'ERROR: ' + T1tofuncMask.cmdline + "\n"
				return 1

			#log that ApplyWarp finished successfully
			print 'SUCCESS: ' +  T1tofuncMask.cmdline + "\n"


	#make func_brain
	ApplyEPIMask=fsl.ImageMaths()
	ApplyEPIMask.inputs.in_file=func
	ApplyEPIMask.inputs.in_file2=EPI_mask
	ApplyEPIMask.inputs.out_file=func_brain
	ApplyEPIMask.inputs.op_string='-mul'

	try:
		ApplyEPIMask.run()
	except:
		print 'ERROR: ' + ApplyEPIMask.cmdline + "\n"
		return 1

	#log that ApplyMask finished successfully
	print 'SUCCESS: ' +  ApplyEPIMask.cmdline + "\n"

	#report successful outcome
	print 'SUCCESS: skullstripEPI\n'
	os.chdir(currentdir)
	outputs=namedtuple("outfiles", ["func_mask","func_brain"])
	return outputs(EPI_mask,func_brain)


#GenTransforms(outDir,ex_func,T1_brain,T1_head,fmap=None,fmapmag=None,fmapmagbrain=None,echospacing=None,pedir=None):
#get rid of fslstandard: use os.environ['FSLDIR'], if this fails exit script
def PreICA(outDir,func,T1_brain,T1_head,fmap=None,fmapmag=None,fmapmagbrain=None,echospacing=None,pedir=None,smooth=True):
	import os
	from collections import namedtuple
	#current directory
	currentdir=os.getcwd()
	if not os.path.isdir(outDir):
		os.makedirs(outDir)
	os.chdir(outDir)
	MCOutputs = MotionCorrect(func=func,outDir=os.path.join(outDir,'mc'))
	GTOutputs = GenTransforms(outDir=os.path.join(outDir,'reg'),ex_func=MCOutputs.mcImgMean,T1_brain=T1_brain,T1_head=T1_head,fmap=fmap,fmapmag=fmapmag,fmapmagbrain=fmapmagbrain,echospacing=echospacing,pedir=pedir)
	SSmcImgOutputs = skullstripEPI(func=MCOutputs.mcImg,T1_brain=T1_brain,T1tofunc_transform=GTOutputs.T1tofunc_transform,EPI_mask=None)
	SSmcImgMeanOutputs = skullstripEPI(func=MCOutputs.mcImgMean,T1_brain=None,T1tofunc_transform=None,EPI_mask=SSmcImgOutputs.func_mask)
	if smooth:
		smoothDir=os.path.join(outDir,'smooth')
		if not os.path.isdir(smoothDir):
			os.makedirs(smoothDir)
		SpatialSmoothOutputs=Spatial_Smooth(func_brain=SSmcImgOutputs.func_brain,func_mask=SSmcImgOutputs.func_mask,func_mean=SSmcImgMeanOutputs.func_brain,outDir=smoothDir)
		PreICAOutsmooth=SpatialSmoothOutputs.smoothed_output
	
	PreICAOutnosmooth=SSmcImgOutputs.func_brain


	os.chdir(currentdir)
	outputs=namedtuple("outfiles", ["func_smooth","func_nosmooth","functoT1_transform","T1toMNI_transform","func_mask","mcImgPar","MNItofunc_warp","functoMNI_warp"])
	return outputs(PreICAOutsmooth,PreICAOutnosmooth,GTOutputs.functoT1_transform,GTOutputs.T1toMNI_transform,SSmcImgOutputs.func_mask,MCOutputs.mcImgPar,GTOutputs.MNItofunc_warp,GTOutputs.functoMNI_warp)

#test: BetaSeries_functions.PreICA('/home/james/RestingState_dev/testout1','/home/james/RestingState_dev/data/pre_sub10_REST_RPI.nii.gz','/home/james/RestingState_dev/data/pre_sub10_MPRAGE-1_RPI_ss.nii.gz','/home/james/RestingState_dev/data/pre_sub10_MPRAGE-1_RPI.nii.gz')

def runICA(fslDir, inFile, outDir, melDirIn, mask, dim, TR):
	""" This function runs MELODIC and merges the mixture modeled thresholded ICs into a single 4D nifti file

	Parameters
	---------------------------------------------------------------------------------
	fslDir:		Full path of the bin-directory of FSL
	inFile:		Full path to the fMRI data file (nii.gz) on which MELODIC should be run
	outDir:		Full path of the output directory
	melDirIn:	Full path of the MELODIC directory in case it has been run before, otherwise define empty string
	mask:		Full path of the mask to be applied during MELODIC
	dim:		Dimensionality of ICA
	TR:		TR (in seconds) of the fMRI data
	
	Output (within the requested output directory)
	---------------------------------------------------------------------------------
	melodic.ica		MELODIC directory
	melodic_IC_thr.nii.gz	merged file containing the mixture modeling thresholded Z-statistical maps located in melodic.ica/stats/ """

	# Import needed modules
	import os
	import commands

	# Define the 'new' MELODIC directory and predefine some associated files
	melDir = os.path.join(outDir,'melodic.ica')
	melIC = os.path.join(melDir,'melodic_IC.nii.gz')
	melICmix = os.path.join(melDir,'melodic_mix')
	melICthr = os.path.join(outDir,'melodic_IC_thr.nii.gz')

	# When a MELODIC directory is specified, check wheter all needed files are present. Otherwise... run MELODIC again
	if (len(melDir) != 0) and os.path.isfile(os.path.join(melDirIn,'melodic_IC.nii.gz')) and os.path.isfile(os.path.join(melDirIn,'melodic_FTmix')) and os.path.isfile(os.path.join(melDirIn,'melodic_mix')):

		print '  - The existing/specified MELODIC directory will be used.'

		# If a 'stats' directory is present (contains thresholded spatial maps) create a symbolic link to the MELODIC directory. Otherwise create specific links and run mixture modeling to obtain thresholded maps.
		if os.path.isdir(os.path.join(melDirIn,'stats')):
			os.symlink(melDirIn,melDir)
		else:
			print '  - The MELODIC directory does not contain the required \'stats\' folder. Mixture modeling on the Z-statistical maps will be run.'
			
			# Create symbolic links to the items in the specified melodic directory
			os.makedirs(melDir)
			for item in os.listdir(melDirIn):
				os.symlink(os.path.join(melDirIn,item),os.path.join(melDir,item))

			# Run mixture modeling
			os.system(' '.join([os.path.join(fslDir,'melodic'),
				'--in=' + melIC,
				'--ICs=' + melIC,
				'--mix=' + melICmix,
				'--outdir=' + melDir,
				'--Ostats --mmthresh=0.5']))
			
	else:
		# If a melodic directory was specified, display that it did not contain all files needed for ICA-AROMA (or that the directory does not exist at all)
		if len(melDirIn) != 0 :
			if not os.path.isdir(melDirIn):
				print '  - The specified MELODIC directory does not exist. MELODIC will be run seperately.'
			else:
				print '  - The specified MELODIC directory does not contain the required files to run ICA-AROMA. MELODIC will be run seperately.'
		
		# Run MELODIC
		os.system(' '.join([os.path.join(fslDir,'melodic'),
			'--in=' + inFile, 
			'--outdir=' + melDir, 
			'--mask=' + mask, 
			'--dim=' + str(dim),
			'--Ostats --nobet --mmthresh=0.5 --report',
			'--tr=' + str(TR)]))

	# Get number of components
	cmd = ' '.join([os.path.join(fslDir,'fslinfo'),
		melIC,
		'| grep dim4 | head -n1 | awk \'{print $2}\''])
	nrICs=int(float(commands.getoutput(cmd)))

	# Merge mixture modeled thresholded spatial maps. Note! In case that mixture modeling did not converge, the file will contain two spatial maps. The latter being the results from a simple null hypothesis test. In that case, this map will have to be used (first one will be empty).
	for i in range(1,nrICs+1):
		# Define thresholded zstat-map file
		zTemp = os.path.join(melDir,'stats','thresh_zstat' + str(i) + '.nii.gz')
		cmd = ' '.join([os.path.join(fslDir,'fslinfo'),
			zTemp,
			'| grep dim4 | head -n1 | awk \'{print $2}\''])
		lenIC=int(float(commands.getoutput(cmd)))

		# Define zeropad for this IC-number and new zstat file
		cmd = ' '.join([os.path.join(fslDir,'zeropad'),
			str(i),
			'4'])
		ICnum=commands.getoutput(cmd)	
		zstat = os.path.join(outDir,'thr_zstat' + ICnum)		

		# Extract last spatial map within the thresh_zstat file
		os.system(' '.join([os.path.join(fslDir,'fslroi'),
			zTemp,		# input
			zstat,		# output
			str(lenIC-1),	# first frame
			'1']))		# number of frames

	# Merge and subsequently remove all mixture modeled Z-maps within the output directory
	os.system(' '.join([os.path.join(fslDir,'fslmerge'),
		'-t',						# concatenate in time
		melICthr,					# output
		os.path.join(outDir,'thr_zstat????.nii.gz')]))	# inputs

	os.system('rm ' + os.path.join(outDir,'thr_zstat????.nii.gz'))

	# Apply the mask to the merged file (in case a melodic-directory was predefined and run with a different mask)
	os.system(' '.join([os.path.join(fslDir,'fslmaths'),
		melICthr,
		'-mas ' + mask,
		melICthr]))

def register2MNI(fslDir, inFile, outFile, affmat, warp):
	""" This function registers an image (or time-series of images) to MNI152 T1 2mm. If no affmat is defined, it only warps (i.e. it assumes that the data has been registerd to the structural scan associated with the warp-file already). If no warp is defined either, it only resamples the data to 2mm isotropic if needed (i.e. it assumes that the data has been registered to a MNI152 template). In case only an affmat file is defined, it assumes that the data has to be linearly registered to MNI152 (i.e. the user has a reason not to use non-linear registration on the data).

	Parameters
	---------------------------------------------------------------------------------
	fslDir:		Full path of the bin-directory of FSL
	inFile:		Full path to the data file (nii.gz) which has to be registerd to MNI152 T1 2mm
	outFile:	Full path of the output file
	affmat:		Full path of the mat file describing the linear registration (if data is still in native space)
	warp:		Full path of the warp file describing the non-linear registration (if data has not been registered to MNI152 space yet)

	Output (within the requested output directory)
	---------------------------------------------------------------------------------
	melodic_IC_mm_MNI2mm.nii.gz	merged file containing the mixture modeling thresholded Z-statistical maps registered to MNI152 2mm """


	# Import needed modules
	import os
	import commands

	# Define the MNI152 T1 2mm template
	fslnobin = fslDir.rsplit('/',2)[0] 
	ref = os.path.join(fslnobin,'data','standard','MNI152_T1_2mm_brain.nii.gz')

	# If the no affmat- or warp-file has been specified, assume that the data is already in MNI152 space. In that case only check if resampling to 2mm is needed
	if (len(affmat) == 0) and (len(warp) == 0):
		# Get 3D voxel size
		pixdim1=float(commands.getoutput('%sfslinfo %s | grep pixdim1 | awk \'{print $2}\'' % (fslDir,inFile) ))
		pixdim2=float(commands.getoutput('%sfslinfo %s | grep pixdim2 | awk \'{print $2}\'' % (fslDir,inFile) ))
		pixdim3=float(commands.getoutput('%sfslinfo %s | grep pixdim3 | awk \'{print $2}\'' % (fslDir,inFile) ))
	
		# If voxel size is not 2mm isotropic, resample the data, otherwise copy the file
		if (pixdim1 != 2) or (pixdim2 != 2) or (pixdim3 !=2 ):
			os.system(' '.join([os.path.join(fslDir,'flirt'),
				' -ref ' + ref,
				' -in ' + inFile,
				' -out ' + outFile,
				' -applyisoxfm 2 -interp trilinear']))
		else:
			os.system('cp ' + inFile + ' ' + outFile)
	
	# If only a warp-file has been specified, assume that the data has already been registered to the structural scan. In that case apply the warping without a affmat
	elif (len(affmat) == 0) and (len(warp) != 0):
		# Apply warp
		os.system(' '.join([os.path.join(fslDir,'applywarp'),
			'--ref=' + ref,
			'--in=' + inFile,
			'--out=' + outFile,
			'--warp=' + warp,
			'--interp=trilinear']))

	# If only a affmat-file has been specified perform affine registration to MNI
	elif (len(affmat) != 0) and (len(warp) == 0):
		os.system(' '.join([os.path.join(fslDir,'flirt'),
			'-ref ' + ref,
			'-in ' + inFile,
			'-out ' + outFile,
			'-applyxfm -init ' + affmat,
			'-interp trilinear']))

	# If both a affmat- and warp-file have been defined, apply the warping accordingly
	else:
		os.system(' '.join([os.path.join(fslDir,'applywarp'),
			'--ref=' + ref,
			'--in=' + inFile,
			'--out=' + outFile,
			'--warp=' + warp,
			'--premat=' + affmat,
			'--interp=trilinear']))


def feature_time_series(melmix, mc):
	""" This function extracts the maximum RP correlation feature scores. It determines the maximum robust correlation of each component time-series with a model of 72 realigment parameters.

	Parameters
	---------------------------------------------------------------------------------
	melmix:		Full path of the melodic_mix text file
	mc:		Full path of the text file containing the realignment parameters
	
	Returns
	---------------------------------------------------------------------------------
	maxRPcorr:	Array of the maximum RP correlation feature scores for the components of the melodic_mix file"""

	# Import required modules
	import numpy as np
	import random

	# Read melodic mix file (IC time-series), subsequently define a set of squared time-series
	mix = np.loadtxt(melmix)
	mixsq = np.power(mix,2)

	# Read motion parameter file
	RP6 = np.loadtxt(mc)

	# Determine the derivatives of the RPs (add zeros at time-point zero)
	RP6_der = np.array(RP6[range(1,RP6.shape[0]),:] - RP6[range(0,RP6.shape[0]-1),:])
	RP6_der = np.concatenate((np.zeros((1,6)),RP6_der),axis=0)

	# Create an RP-model including the RPs and its derivatives
	RP12 = np.concatenate((RP6,RP6_der),axis=1)

	# Add the squared RP-terms to the model
	RP24 = np.concatenate((RP12,np.power(RP12,2)),axis=1)

	# Derive shifted versions of the RP_model (1 frame for and backwards)
	RP24_1fw = np.concatenate((np.zeros((1,24)),np.array(RP24[range(0,RP24.shape[0]-1),:])),axis=0)
	RP24_1bw = np.concatenate((np.array(RP24[range(1,RP24.shape[0]),:]),np.zeros((1,24))),axis=0)

	# Combine the original and shifted mot_pars into a single model
	RP_model = np.concatenate((RP24,RP24_1fw,RP24_1bw),axis=1)

	# Define the column indices of respectively the squared or non-squared terms
	idx_nonsq = np.array(np.concatenate((range(0,12), range(24,36), range(48,60)),axis=0))
	idx_sq = np.array(np.concatenate((range(12,24), range(36,48), range(60,72)),axis=0))

	# Determine the maximum correlation between RPs and IC time-series
	nSplits=int(1000)
	maxTC = np.zeros((nSplits,mix.shape[1]))
	for i in range(0,nSplits):
		# Get a random set of 90% of the dataset and get associated RP model and IC time-series matrices
		idx = np.array(random.sample(range(0,mix.shape[0]),int(round(0.9*mix.shape[0]))))
		RP_model_temp = RP_model[idx,:]
		mix_temp = mix[idx,:]
		mixsq_temp = mixsq[idx,:]

		# Calculate correlation between non-squared RP/IC time-series
		RP_model_nonsq = RP_model_temp[:,idx_nonsq]
		cor_nonsq = np.array(np.zeros((mix_temp.shape[1],RP_model_nonsq.shape[1])))
		for j in range(0,mix_temp.shape[1]):
			for k in range(0,RP_model_nonsq.shape[1]):
				cor_temp = np.corrcoef(mix_temp[:,j],RP_model_nonsq[:,k])
				cor_nonsq[j,k] = cor_temp[0,1]

		# Calculate correlation between squared RP/IC time-series
		RP_model_sq = RP_model_temp[:,idx_sq]
		cor_sq = np.array(np.zeros((mix_temp.shape[1],RP_model_sq.shape[1])))
		for j in range(0,mixsq_temp.shape[1]):
			for k in range(0,RP_model_sq.shape[1]):
				cor_temp = np.corrcoef(mixsq_temp[:,j],RP_model_sq[:,k])
				cor_sq[j,k] = cor_temp[0,1]

		# Combine the squared an non-squared correlation matrices
		corMatrix = np.concatenate((cor_sq,cor_nonsq),axis=1)

		# Get maximum absolute temporal correlation for every IC
		corMatrixAbs = np.abs(corMatrix)
		maxTC[i,:] = corMatrixAbs.max(axis=1)

	# Get the mean maximum correlation over all random splits
	maxRPcorr = maxTC.mean(axis=0)

	# Return the feature score
	return maxRPcorr

def feature_frequency(melFTmix, TR):
	""" This function extracts the high-frequency content feature scores. It determines the frequency, as fraction of the Nyquist frequency, at which the higher and lower frequencies explain half of the total power between 0.01Hz and Nyquist. 
	
	Parameters
	---------------------------------------------------------------------------------
	melFTmix:	Full path of the melodic_FTmix text file
	TR:		TR (in seconds) of the fMRI data (float)
	
	Returns
	---------------------------------------------------------------------------------
	HFC:		Array of the HFC ('High-frequency content') feature scores for the components of the melodic_FTmix file"""

	# Import required modules
	import numpy as np

	# Determine sample frequency
	Fs = 1/TR

	# Determine Nyquist-frequency
	Ny = Fs/2
		
	# Load melodic_FTmix file
	FT=np.loadtxt(melFTmix)


	# Determine which frequencies are associated with every row in the melodic_FTmix file  (assuming the rows range from 0Hz to Nyquist)
	f = Ny*(np.array(range(1,FT.shape[0]+1)))/(FT.shape[0])


	# Only include frequencies higher than 0.01Hz
	fincl = np.squeeze(np.array(np.where( f > 0.01 )))
	FT=FT[fincl,:]
	f=f[fincl]

	# Set frequency range to [0-1]
	f_norm = (f-0.01)/(Ny-0.01)

	# For every IC; get the cumulative sum as a fraction of the total sum
	fcumsum_fract = np.cumsum(FT,axis=0)/ np.sum(FT,axis=0)

	# Determine the index of the frequency with the fractional cumulative sum closest to 0.5
	idx_cutoff=np.argmin(np.abs(fcumsum_fract-0.5),axis=0)

	# Now get the fractions associated with those indices index, these are the final feature scores
	HFC = f_norm[idx_cutoff]
		 
	# Return feature score
	return HFC

def feature_spatial(fslDir, tempDir, aromaDir, melIC):
	""" This function extracts the spatial feature scores. For each IC it determines the fraction of the mixture modeled thresholded Z-maps respecitvely located within the CSF or at the brain edges, using predefined standardized masks.

	Parameters
	---------------------------------------------------------------------------------
	fslDir:		Full path of the bin-directory of FSL
	tempDir:	Full path of a directory where temporary files can be stored (called 'temp_IC.nii.gz')
	aromaDir:	Full path of the ICA-AROMA directory, containing the mask-files (mask_edge.nii.gz, mask_csf.nii.gz & mask_out.nii.gz) 
	melIC:		Full path of the nii.gz file containing mixture-modeled threholded (p>0.5) Z-maps, registered to the MNI152 2mm template
	
	Returns
	---------------------------------------------------------------------------------
	edgeFract:	Array of the edge fraction feature scores for the components of the melIC file
	csfFract:	Array of the CSF fraction feature scores for the components of the melIC file"""
	
	# Import required modules
	import numpy as np
	import os
	import commands

	currentdir=os.getcwd()
	os.chdir(aromaDir)
	# Get the number of ICs
	numICs = int(commands.getoutput('%sfslinfo %s | grep dim4 | head -n1 | awk \'{print $2}\'' % (fslDir, melIC) ))

	# Loop over ICs
	edgeFract=np.zeros(numICs)
	csfFract=np.zeros(numICs)
	for i in range(0,numICs):
		# Define temporary IC-file
		tempIC = os.path.join(tempDir,'temp_IC.nii.gz')

		# Extract IC from the merged melodic_IC_thr2MNI2mm file
		os.system(' '.join([os.path.join(fslDir,'fslroi'),
			melIC,
			tempIC,
			str(i),
			'1']))

		# Change to absolute Z-values
		os.system(' '.join([os.path.join(fslDir,'fslmaths'),
			tempIC,
			'-abs',
			tempIC]))
		
		# Get sum of Z-values within the total Z-map (calculate via the mean and number of non-zero voxels)
		totVox = int(commands.getoutput(' '.join([os.path.join(fslDir,'fslstats'),
							tempIC,
							'-V | awk \'{print $1}\''])))
		
		if not (totVox == 0):
			totMean = float(commands.getoutput(' '.join([os.path.join(fslDir,'fslstats'),
							tempIC,
							'-M'])))
		else:
			print '     - The spatial map of component ' + str(i+1) + ' is empty. Please check!'
			totMean = 0

		totSum = totMean * totVox
		
		# Get sum of Z-values of the voxels located within the CSF (calculate via the mean and number of non-zero voxels)
		csfVox = int(commands.getoutput(' '.join([os.path.join(fslDir,'fslstats'),
							tempIC,
							'-k mask_csf.nii.gz',
							'-V | awk \'{print $1}\''])))

		if not (csfVox == 0):
			csfMean = float(commands.getoutput(' '.join([os.path.join(fslDir,'fslstats'),
							tempIC,
							'-k mask_csf.nii.gz',
							'-M'])))
		else:
			csfMean = 0

		csfSum = csfMean * csfVox	

		# Get sum of Z-values of the voxels located within the Edge (calculate via the mean and number of non-zero voxels)
		edgeVox = int(commands.getoutput(' '.join([os.path.join(fslDir,'fslstats'),
							tempIC,
							'-k mask_edge.nii.gz',
							'-V | awk \'{print $1}\''])))
		if not (edgeVox == 0):
			edgeMean = float(commands.getoutput(' '.join([os.path.join(fslDir,'fslstats'),
							tempIC,
							'-k mask_edge.nii.gz',
							'-M'])))
		else:
			edgeMean = 0
		
		edgeSum = edgeMean * edgeVox

		# Get sum of Z-values of the voxels located outside the brain (calculate via the mean and number of non-zero voxels)
		outVox = int(commands.getoutput(' '.join([os.path.join(fslDir,'fslstats'),
							tempIC,
							'-k mask_out.nii.gz',
							'-V | awk \'{print $1}\''])))
		if not (outVox == 0):
			outMean = float(commands.getoutput(' '.join([os.path.join(fslDir,'fslstats'),
							tempIC,
							'-k mask_out.nii.gz',
							'-M'])))
		else:
			outMean = 0
		
		outSum = outMean * outVox

		# Determine edge and CSF fraction
		if not (totSum == 0):
			edgeFract[i] = (outSum + edgeSum)/(totSum - csfSum)
			csfFract[i] = csfSum / totSum
		else:
			edgeFract[i]=0
			csfFract[i]=0

	# Remove the temporary IC-file
	os.remove(tempIC)
	os.chdir(currentdir)
	# Return feature scores
	return edgeFract, csfFract

def classification(outDir, maxRPcorr, edgeFract, HFC, csfFract):
	""" This function classifies a set of components into motion and non-motion components based on four features; maximum RP correlation, high-frequency content, edge-fraction and CSF-fraction

	Parameters
	---------------------------------------------------------------------------------
	outDir:		Full path of the output directory
	maxRPcorr:	Array of the 'maximum RP correlation' feature scores of the components
	edgeFract:	Array of the 'edge fraction' feature scores of the components
	HFC:		Array of the 'high-frequency content' feature scores of the components
	csfFract:	Array of the 'CSF fraction' feature scores of the components

	Return
	---------------------------------------------------------------------------------
	motionICs	Array containing the indices of the components identified as motion components

	Output (within the requested output directory)
	---------------------------------------------------------------------------------
	classified_motion_ICs.txt	A text file containing the indices of the components identified as motion components """

	# Import required modules
	import numpy as np
	import os
	import commands

	# Classify the ICs as motion or non-motion

	# Define criteria needed for classification (thresholds and hyperplane-parameters)
	thr_csf = 0.10
	thr_HFC = 0.35
	hyp = np.array([-19.9751070082159, 9.95127547670627, 24.8333160239175])
	
	# Project edge & maxRPcorr feature scores to new 1D space
	x = np.array([maxRPcorr, edgeFract])
	proj = hyp[0] + np.dot(x.T,hyp[1:])

	# Classify the ICs
	motionICs = np.squeeze(np.array(np.where((proj > 0) + (csfFract > thr_csf) + (HFC > thr_HFC))))

	# Put the feature scores in a text file
	np.savetxt(os.path.join(outDir,'feature_scores.txt'),np.vstack((maxRPcorr,edgeFract,HFC,csfFract)).T)

	# Put the indices of motion-classified ICs in a text file
	txt = open(os.path.join(outDir,'classified_motion_ICs.txt'),'w')
	if len(motionICs) != 0:
		txt.write(','.join(['%.0f' % num for num in (motionICs+1)]))
	txt.close()

	# Create a summary overview of the classification
	txt = open(os.path.join(outDir,'classification_overview.txt'),'w')
	txt.write('IC' + '\t' +  'Motion/noise' + '\t' +  'maximum RP correlation' + '\t' +  'Edge-fraction' + '\t\t' +  'High-frequency content' + '\t' + 'CSF-fraction')
	txt.write('\n')
	for i in range(0,len(csfFract)):
		if (proj[i] > 0) or (csfFract[i] > thr_csf) or (HFC[i] > thr_HFC):
			classif="True"
		else:
			classif="False"
		txt.write('%.0f\t%s\t\t%.2f\t\t\t%.2f\t\t\t%.2f\t\t\t%.2f\n' % (i+1, classif, maxRPcorr[i], edgeFract[i], HFC[i], csfFract[i]))
	txt.close()

	return motionICs

def denoising(fslDir, inFile, outDir, melmix, denType, denIdx):
	""" This function classifies the ICs based on the four features; maximum RP correlation, high-frequency content, edge-fraction and CSF-fraction

	Parameters
	---------------------------------------------------------------------------------
	fslDir:		Full path of the bin-directory of FSL
	inFile:		Full path to the data file (nii.gz) which has to be denoised
	outDir:		Full path of the output directory
	melmix:		Full path of the melodic_mix text file
	denType:	Type of requested denoising ('aggr': aggressive, 'nonaggr': non-aggressive, 'both': both aggressive and non-aggressive 
	denIdx:		Indices of the components that should be regressed out

	Output (within the requested output directory)
	---------------------------------------------------------------------------------
	denoised_func_data_<denType>.nii.gz:		A nii.gz file of the denoised fMRI data"""

	# Import required modules
	import os
	import numpy as np

	# Check if denoising is needed (i.e. are there components classified as motion)
	check = len(denIdx) > 0

	if check==1:
		# Put IC indices into a char array
		denIdxStr = np.char.mod('%i',(denIdx+1))

		# Non-aggressive denoising of the data using fsl_regfilt (partial regression), if requested
		if (denType == 'nonaggr') or (denType == 'both'):		
			os.system(' '.join([os.path.join(fslDir,'fsl_regfilt'),
				'--in=' + inFile,
				'--design=' + melmix,
				'--filter="' + ','.join(denIdxStr) + '"',
				'--out=' + os.path.join(outDir,'denoised_func_data_nonaggr.nii.gz')]))

		# Aggressive denoising of the data using fsl_regfilt (full regression)
		if (denType == 'aggr') or (denType == 'both'):
			os.system(' '.join([os.path.join(fslDir,'fsl_regfilt'),
				'--in=' + inFile,
				'--design=' + melmix,
				'--filter="' + ','.join(denIdxStr) + '"',
				'--out=' + os.path.join(outDir,'denoised_func_data_aggr.nii.gz'),
				'-a']))
	else:
		print "  - None of the components was classified as motion, so no denoising is applied (a symbolic link to the input file will be created)."
		if (denType == 'nonaggr') or (denType == 'both'):
			os.symlink(inFile,os.path.join(outDir,'denoised_func_data_nonaggr.nii.gz'))
		if (denType == 'aggr') or (denType == 'both'):
			os.symlink(inFile,os.path.join(outDir,'denoised_func_data_aggr.nii.gz'))

#import BetaSeries_functions
#test:  BetaSeries_functions.GenTransforms('/usr/share/fsl/5.0/','/home/james/RestingState_dev/','/home/james/RestingState_dev/data/pre_sub10_REST_RPI.nii.gz','/home/james/RestingState_dev/data/pre_sub10_MPRAGE-1_RPI_ss.nii.gz','/home/james/RestingState_dev/data/pre_sub10_MPRAGE-1_RPI.nii.gz')

# BetaSeries_functions.GenTransforms('/usr/share/fsl/data/standard/','/home/james/RestingState_dev/reg','/home/james/RestingState_dev/test/mc/mcImgmean.nii.gz','/home/james/RestingState_dev/data/pre_sub10_MPRAGE-1_RPI_ss.nii.gz','/home/james/RestingState_dev/data/pre_sub10_MPRAGE-1_RPI.nii.gz')


##Post ICA
def TemporalFilter(denoised_func,outDir):
	import os
	import string
	import nipype.interfaces.fsl as fsl
	currentdir=os.getcwd()
	
	#go to outdir
	if not os.path.isdir(outDir):
		os.makedirs(outDir)
	os.chdir(outDir)
	Tempfilt_name=os.path.basename(denoised_func).replace('.nii.gz','_temp_filt.nii.gz')
	Tempfilt_output=os.path.join(outDir,Tempfilt_name)

	#highpass calculated from TR=2.0 seconds
	TemporalFilterCmd=fsl.TemporalFilter(in_file=denoised_func,highpass_sigma=25.0,out_file=Tempfilt_output)
	TemporalFilterCmd.run()
	os.chdir(currentdir)
	return Tempfilt_output

def NuisanceRegression(filtered_func,Nrois,MNItofuncWarp,outdir,motion=False,eig=False):
	import os
	import string
	import nipype.interfaces.fsl as fsl
	import numpy as np
	import nibabel as nib
	import subprocess

	#get the current directory:
	currentdir=os.getcwd()
	#go the outdir to make sure "intermediary" files are dropped off in the same place as "final" files
	os.chdir(outdir)
	#read the nuisance rois txt file
	with open(Nrois,'r') as txt:
		nrois=txt.read().strip().split('\n')
	#initialize overall numpy array
	func=nib.load(filtered_func)
	length_ts=func.shape[3]
	num_nrois=len(nrois)
	#np.empty([num_nrois,length_ts])
	if not eig:
		nrois_norm_ts_tot=os.path.join(outdir,'nrois_norm_ts.txt')
		nrois_norm_ts_mat=os.path.join(outdir,'nrois_norm_ts.mat')
	if eig:
		nrois_norm_ts_tot=os.path.join(outdir,'_eig_nrois_norm_ts.txt')
		nrois_norm_ts_mat=os.path.join(outdir,'_eig_nrois_norm_ts.mat')

	NuisanceReg_func=os.path.join(outdir,'NuisanceReg.nii.gz')
	for index,nroi in enumerate(nrois,start=0):
		basenroi=os.path.basename(nroi)
		subnroi=os.path.join(outdir,basenroi)
		nroi_ts=string.replace(subnroi,'.nii.gz','_ts.txt')
		nroi_bin=string.replace(subnroi,'.nii.gz','_bin.nii.gz')
		nroi_norm_ts=string.replace(subnroi,'.nii.gz','_norm_ts.txt')
		#nrois_norm_ts_tot=string.replace(subnroi,'.nii.gz','_norm_ts_tot.txt')
		#nrois_norm_ts_mat=string.replace(subnroi,'.nii.gz','_norm_ts_tot.mat')
		#move nroi to subject space
		stdnroi2subnroi=fsl.ApplyWarp(in_file=nroi,out_file=subnroi,ref_file=filtered_func,field_file=MNItofuncWarp)
		stdnroi2subnroi.run()

		#binarize mask
		nroiBinCmd=fsl.ImageMaths(in_file=subnroi,out_file=nroi_bin,op_string='-thr 0.5 -bin')
		nroiBinCmd.run()
		
		#extract the time series from the mask
		if not eig:
			ts_extraction=fsl.ImageMeants(in_file=filtered_func,mask=nroi_bin,out_file=nroi_ts)
			ts_extraction.run()
		elif eig:
			nroi_base=string.replace(subnroi,'.nii.gz','')
			nroi_ts=string.replace(subnroi,'.nii.gz','00.1D')
			roi_bin=string.replace(subnroi,'.nii.gz','_eig_bin.nii.gz')
			nroi_norm_ts=string.replace(subnroi,'.nii.gz','_eig_norm_ts.txt')
			subprocess.check_output("3dpc -prefix %s -pcsave 1 -nscal -mask %s %s" % (subnroi_base,nroi_bin,filtered_func),shell=True)
		else:
			print 'Error! eig not set'
			return 1

		#read in the ts file
		ts_arr=np.loadtxt(nroi_ts)

		#norm the array
		ts_norm=(ts_arr - ts_arr.mean()) / ts_arr.std()
		
		#write normed array to file
		np.savetxt(nroi_norm_ts,ts_norm,'%3.10f')

		#save in larger array
		if index == 0:
			ts_norm_tot=ts_norm
		else:
			ts_norm_tot=np.column_stack((ts_norm_tot,ts_norm))

	#print out the total array to a file
	np.savetxt(nrois_norm_ts_tot,ts_norm_tot,'%3.10f')

	#transfer the text file to be readable by fsl
	os.system('Text2Vest %s %s' % (nrois_norm_ts_tot,nrois_norm_ts_mat))

	#return the normed list (for use to apply the regres,sion)
	normed_list=[string.replace(nroi,'.nii.gz','_norm_ts.txt') for nroi in nrois]

	#perform nuisance regression
	nuisanceGLM=fsl.GLM(in_file=filtered_func,design=nrois_norm_ts_mat,out_res_name=NuisanceReg_func)
	nuisanceGLM.run()

	#go back to the original directory
	os.chdir(currentdir)

	return NuisanceReg_func

def BetaSeries(outdir,whichEVs,numrealev,tempderiv=None):
	#assumes the requisite files exist in an fsl "friendly" way
	#use symlinks in the main script to connect the images
	#i.e. mask.nii.gz, filtered_func_data.nii.gz, and custom_timing_files exist
	import os
	import string
	import nipype.interfaces.fsl as fsl
	import subprocess
	#tmp addition directory of pybetaseries
	subprocess.check_output("pybetaseries.py --fsldir %s --whichevs %s --numrealev %s %s -LSS" % (outdir,str(whichEVs).strip('[]').replace(',',''),numrealev,tempderiv),shell=True)
	print "BetaSeries Successful!"
	
def MakeModel(func,EVs,outDir):
	#notes: need to make the outdir have design.mat and folder "custom_timing_files with properly named evs inside"
	#so clean up the evs in the main script before calling this function
	import os
	import nipype.algorithms.modelgen as modelgen
	import nipype.interfaces.fsl.model as model
	import shutil
	#ugg: I dislike changing directories for the sake of output!!!
	currentdir=os.getcwd()
	#outDir=os.path.dirname(func)
	design_fsf=os.path.join(outDir,'design.fsf')
	design_mat=os.path.join(outDir,'design.mat')

	EVDir=os.path.join(outDir,'custom_timing_files')
	if not os.path.isdir(EVDir):
		os.makedirs(EVDir)
	for index,EV in enumerate(EVs,start=1):
		EVoutfile=os.path.join(EVDir,'ev%s.txt' % (index))
		shutil.copyfile(EV,EVoutfile)


	os.chdir(outDir)
	ModelSpec=modelgen.SpecifyModel()
	ModelSpec.inputs.event_files=EVs #a list of 3-column txt files
	ModelSpec.inputs.functional_runs=func #the preprocessed functional volume
	ModelSpec.inputs.high_pass_filter_cutoff=100 #make this user changable?
	ModelSpec.inputs.input_units='secs'
	ModelSpec.inputs.time_repetition=2 #change this to read from the data?
	subinfo=ModelSpec.run()

	FirstLevel=model.Level1Design()
	FirstLevel.inputs.interscan_interval=2
	FirstLevel.inputs.bases={'dgamma':{'derivs':True}}
	FirstLevel.inputs.session_info=subinfo.outputs.session_info
	FirstLevel.inputs.model_serial_correlations=True
	fsf_res=FirstLevel.run()

	FEAT=model.FEATModel()
	FEAT.inputs.ev_files=EVs
	FEAT.inputs.fsf_file=fsf_res.outputs.fsf_files
	design_res=FEAT.run() #get the design file for pybetaseries

	#rename the files for pybetaseries
	os.rename(fsf_res.outputs.fsf_files,design_fsf)
	os.rename(design_res.outputs.design_file,design_mat)

	os.chdir(currentdir)



def Normalize(Image):
	import os
	import nipype.interfaces.fsl as fsl
	currentdir=os.getcwd()
	outDir=os.path.dirname(Image)
	os.chdir(outDir)


	#outfile names
	ImageMean=Image.replace('.nii.gz','_mean.nii.gz')
	ImageSD=Image.replace('.nii.gz','_sd.nii.gz')
	ImageNorm=Image.replace('.nii.gz','_norm.nii.gz')

	#get the standard deviation
	SDimg=fsl.maths.StdImage(in_file=Image,out_file=ImageSD)
	SDimg.run()

	#get the mean
	MEANimg=fsl.MeanImage(in_file=Image,out_file=ImageMean)
	MEANimg.run()

	#Normalize the data to have a mean of 0 and an SD of 1.
	NORMimg=fsl.maths.MultiImageMaths(in_file=Image,op_string='-sub %s -div %s',operand_files=[ImageMean,ImageSD],out_file=ImageNorm)
	NORMimg.run()

	return ImageNorm


def SeedCorrelate(EVLSS,seed,Seed_Outdir,MNItofuncWarp,functoMNIwarp,eig=False):
	import os
	import nipype.interfaces.fsl as fsl
	import nipype.interfaces.afni as afni
	import numpy as np
	import string
	if eig:
		import subprocess

	currentdir=os.getcwd()
	if not os.path.isdir(Seed_Outdir):
		os.makedirs(Seed_Outdir)
	os.chdir(Seed_Outdir)

	#outfile names
	baseseed=os.path.basename(seed)
	subseed=os.path.join(Seed_Outdir,baseseed)
	subseed_bin=string.replace(subseed,'.nii.gz','_bin.nii.gz')
	seed_ts=string.replace(subseed,'.nii.gz','_ts.txt')
	seed_norm_ts=string.replace(subseed,'.nii.gz','_norm_ts.txt')
	seed_rscore=string.replace(subseed,'.nii.gz','_rscore.nii.gz')
	seed_rscore_MNI=string.replace(subseed,'.nii.gz','_rscore_MNI.nii.gz')
	seed_zscore_MNI=string.replace(subseed,'.nii.gz','_zscore_MNI.nii.gz')
	#Normalize the dataset:
	EVLSS_norm = Normalize(EVLSS)

	#move nroi to subject space
	stdseed2subseed=fsl.ApplyWarp(in_file=seed,out_file=subseed,ref_file=EVLSS,field_file=MNItofuncWarp)
	stdseed2subseed.run()

	#binarize mask
	seedBinCmd=fsl.ImageMaths(in_file=subseed,out_file=subseed_bin,op_string='-thr 0.5 -bin')
	seedBinCmd.run()
	
	#extract the time series from the mask
	if not eig:
		ts_extraction=fsl.ImageMeants(in_file=EVLSS_norm,mask=subseed_bin,out_file=seed_ts)
		ts_extraction.run()
	elif eig:
		 seedtsbase=string.replace(subseed,'.nii.gz','')
		 seed_ts=string.replace(subseed,'.nii.gz','00.1D')
		 seed_norm_ts=string.replace(subseed,'.nii.gz','_norm_eig_ts.txt')
		 seed_rscore=string.replace(subseed,'.nii.gz','_eig_rscore.nii.gz')
		 seed_rscore_MNI=string.replace(subseed,'.nii.gz','_eig_rscore_MNI.nii.gz')
		 seed_zscore_MNI=string.replace(subseed,'.nii.gz','_eig_zscore_MNI.nii.gz')
		 subprocess.check_output("3dpc -prefix %s -pcsave 1 -nscal -mask %s %s" % (seedtsbase,subseed_bin,EVLSS_norm),shell=True)
	else:
		print 'eig not define: error!'
		return 1

	#read in the ts file
	ts_arr=np.loadtxt(seed_ts)

	#norm the array
	ts_norm=(ts_arr - ts_arr.mean()) / ts_arr.std()
	
	#write normed array to file
	np.savetxt(seed_norm_ts,ts_norm,'%3.10f')

	#correlate the normed ts with the normed image
	correlate=afni.TCorr1D(xset=EVLSS_norm,y_1d=seed_norm_ts,out_file=seed_rscore)
	correlate.run()

	#register the image to MNI space
	sub2MNI=fsl.ApplyWarp(in_file=seed_rscore,out_file=seed_rscore_MNI,ref_file=fsl.Info.standard_image('MNI152_T1_2mm_brain.nii.gz'),field_file=functoMNIwarp)
	sub2MNI.run()

	#transfer the rscore to a zscore
	r2z_transform=afni.Calc(in_file_a=seed_rscore_MNI,expr='log((1+a)/(1-a))/2',out_file=seed_zscore_MNI)
	r2z_transform.run()

	#change directory back to original
	os.chdir(currentdir)
	#this is it?


#Notes:
#Don't do prewhitening for a seed based analysis?
#https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1203&L=fsl&D=0&P=300924