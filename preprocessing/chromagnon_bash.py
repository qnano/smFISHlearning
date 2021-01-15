##################################################################
# Chromagnon bash
#
# -Runs chromagnon channel alignment tool on an entire directory
#
# Mode 1- each targ file has a ref file
# -requires ref and targ directory
# -requires ref and targ images to have identical names
#
# Mode 2- one ref file and a directory of targ files
# Requires .csv ref file and targ directory
#
# --USAGE--
# - choose mode by calling either function
# - python ~/src/Python_Scripts/chromagnon_bash.py
#
##################################################################


import subprocess
import os

def ref_file_targ_file():

	# specify directory of reference and target images
	indir = '/usr/people/bioc1301/data/AdultBrain_smFISH_MASTER/250220_280220_smFISH_MBONs/20200331_repeat'
	ref_dir = os.path.join(indir, 'cal')
	ref_files = os.listdir(ref_dir)
	targ_dir = os.path.join(indir, 'image')
	targ_files = os.listdir(targ_dir)

	for file in ref_files:
            if file.endswith('.tif'):
    	        targ = os.path.join(targ_dir, file)
	        ref = os.path.join(ref_dir, file)

	        # call chromagnon from bash
	        export_command = ["chromagnon", targ, "-R", ref, "-E", "dv"]
                print("calling ", export_command)
                try:
	            subprocess.check_call(export_command)
                except CalledProcessError as e:
                    print('failed to run %s: %e' % (export_command, str(e)))

def ref_file_targ_file_rename():

	# specify directory of reference and target images
	ref_dir = ('/usr/people/bioc1301/data/20191018_R66C08_smFISH/cal/complete')
	ref_files = os.listdir(ref_dir)
	targ_dir = ('/usr/people/bioc1301/data/20191018_R66C08_smFISH/image')
	targ_files = os.listdir(targ_dir)

	for file in ref_files:
		if file.endswith('.csv'):
			print ('aligning:',file)
			ref = os.path.join(ref_dir, file)
			targ = os.path.join(targ_dir, file[:-15])
			# call chromagnon from bash
			export_command = ["chromagnon", targ, "-R", ref, "-E", "dv"]
			#export_command = "chromagnon "+targ+" -R "+ref+" -E dv"
			#popen = subprocess.Popen(export_command, shell=True, stdout=subprocess.PIPE)
			#out, err = popen.communicate()
			#output = out.split(':')[0]
			#print (out)
			#print("calling ", export_command)
                	print("calling ", export_command)
			try:
                    	    subprocess.check_call(export_command)
                	except CalledProcessError as e:
                    	    print('failed to run %s: %e' % (export_command, str(e)))


def one_ref_file():

	# specify directory of images
	targ_dir = ('/usr/people/bioc1301/data/Olympus_FV3000/20190621_eIF4eGFP_Kstim_smFISH/oir/')
	targ_files = os.listdir(targ_dir)

	for file in targ_files:
		if file.endswith('.tif'):
    			print ('aligning:', file)
    			export_command = "chromagnon "+file+" -R /usr/people/bioc1301/data/Olympus_FV3000/20190621_eIF4eGFP_Kstim_smFISH/20190621_eIF4eGFP_msp670_syp568_HRP_viol_cal.tif.chromagnon.csv -E dv"
    			popen = subprocess.Popen(export_command, shell=True, stdout=subprocess.PIPE)
    			out, err = popen.communicate()
    			output = str(out).split(':')[0]
    			print (out)


#ref_file_targ_file_rename()
ref_file_targ_file()
#one_ref_file()
