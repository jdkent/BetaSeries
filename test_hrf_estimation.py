#!/usr/bin/env python

#required packages for setup/main script
import argparse

def hrf_estimation(tstxt,TR,stimtimes,pretime,posttime):
	"""################
	   #function: hrf_estimation
	   ################
	   #purpose: estimates a hemodynamic response from a txt file timeseries
	   ################
	   #preconditions: afni is installed
	   ################
	   #input: tstxt,txt, the 1D file containing the timeseries you want estimated
	   #	   TR,integer, the repeition time between volumes
	   #	   stimtimes,list, a list containing the text files that represent the stimulus onset timeseries
	   #	   pretime,integer, the number of seconds before the stimulus you want estimated
	   #	   posttime,integer, the number of seconds after the stimulus you want estimated 
	   ################
	   #output: *_resp.1D for each condition, gives estimate of HRF shape at each TR (not each second)
	   #		***other output not documented (subject to change)
	   ################"""
	import subprocess
	import os

	#creating variables to pass into AFNI
	numstimts=len(stimtimes)
	numpoints=((posttime-pretime)/TR)+1 # how many timepoints in the hrf are being estimated?
	tstxtname=os.path.basename(os.path.splitext(tstxt)[0])
	condnames=[ os.path.basename(os.path.splitext(x)[0]) for x in stimtimes ]
	stimts_input=[] # initialize list (I don't think this is very pythonic)
	for num,txt in enumerate(zip(stimtimes,condnames)):

		stimfile=open(txt[0],"r")
		stimlist=[ float(x) for x in stimfile.read().split('\n') if x ]
		if len(stimlist) is 0 or len(stimlist) is 1 and stimlist[0] == 0:
			print "stim file either is empty or only has a zero: " + txt[0]
			numstimts-=1
		else:
			stimts_input.append("-stim_times %s %s 'TENT(%s,%s,%s)' -stim_label %s %s -iresp %s %s_resp" % (num+1, txt[0], pretime, posttime, numpoints, num+1, txt[1], num+1, txt[1]))

	#treat the list as a long string
	stimts_input_string=' '.join(str(x) for x in stimts_input)

	#condtruct the entire afni command
	#caution, I'm including the GOFORIT 3 command becuase some the error files only have 1 value, and this makes it very hard to estimate the HRF using FIR with only one observation.
	#Since I don't care about error trials, I'm telling 3dDeconvolve to run anyways.
	afni_cmd="3dDeconvolve -GOFORIT 3 -input1D %s -polort A -TR_1D %s -num_stimts %s %s -tout -fout -bucket %s_bucket" % (tstxt,TR,numstimts,stimts_input_string,tstxtname)

	#print what we are running
	print "Running:\n" + afni_cmd + "\n"

	#actually run the afni command
	subprocess.check_output(afni_cmd,shell=True)
#########################################################################




#set up arguments to be read in from the command line
parser = argparse.ArgumentParser(description='wrapper to use AFNI\'s 3dDeconvolve to estimate the hrf from a 1d txt file')

parser.add_argument('-t','--tstxt',dest="tstxt",required=True, help='the timeseries file listed in 1D format (e.g. one number per line)')
parser.add_argument('-T','--TR',dest="TR",required=True, help='The repitition time (e.g. time to collect a volume), note: currently only supports integers)')
parser.add_argument('-s','--stimfile',dest="stimfiles",required=True,nargs='+', help='One column stimulus onset files, one number per line, and separate multiple text files with a space (e.g. -s stim1.txt stim2.txt ... )')
parser.add_argument('-p','--pretime',dest="pretime",required=True, help='The number of seconds before the stimulus you are interested in plotting (note: final output only gives datapoints at the resolution of TR)')
parser.add_argument('-P','--postime',dest="posttime",required=True, help='The number of seconds after the stimulus you are interested in plotting (note: final output only gives datapoints at the resolution of TR)')

#parse the arguments
args = parser.parse_args()

#main script
hrf_estimation(tstxt=args.tstxt,TR=int(args.TR),stimtimes=args.stimfiles,pretime=int(args.pretime),posttime=int(args.posttime))

#example call
#test_hrf_estimation.py -t VTA_bin00.1D -T 2 -s s10_con_1column.txt s10_inc_1column.txt s10_neu_1column.txt s10_errors_1column.txt -p -4 -P 20

#type test_hrf_estimation.py -h or --help for help