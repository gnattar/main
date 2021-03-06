#!/usr/bin/env python

####################################################################################
#
# S Peron 2011 Sept
#
#   1) looks for ready_to_run files, meaking breakup is done and time for meat
#   2) kicks off a bunch of cluster processes 
#
####################################################################################
import sys, struct, os, re, string, shutil, random;
import xml.dom.minidom, time;

# definitely (?) fixed -- where you go huntin for dem whiska files
root_search_path = "/groups/magee/mageelab/GR_dm11/whisker/reprocessed";
ready_filename = "ready_for_linking";
qsub_path = "/sge/8.0.1p4/bin/lx-amd64/qsub";
par_execute_path = "/groups/magee/mageelab/GR_dm11/imreg_on_cluster_GR/anal/par/par_execute_cluster";
logs_path = "/groups/magee/mageelab/GR_dm11/whisker/reprocessed/logs";
n_jobs_per_pardir = 100;

# path passed?
if (len(sys.argv) > 1):
  root_search_path = sys.argv[1];

#
# This will take a single drectory of parfiles and generate many jobs for it
#
def process_pardir(parpath):
	uid = "%08d" % (random.random()*100000000) # so we don't collide
	running_file_full = parpath + "/" + ready_filename;

	# qsub 
	for f in range(1,n_jobs_per_pardir):
		fs = "%d" % (f);
		jobname = "swtlink-" + fs;
		dispatch_cmd = qsub_path + ' -N ' + jobname + ' -pe batch 1 -j y -o /dev/null -b y -cwd -V "' + par_execute_path + ' ' + parpath + ' > ' + logs_path  + '/swtlink-' + uid + '-' + fs + '.log"';
		print dispatch_cmd;
		os.system(dispatch_cmd);
		time.sleep(5); # wait 5 seconds before kickoff of next ...

	# Remove ready file
	shutil.copyfile(running_file_full,running_file_full + ".done");
	os.remove(running_file_full);


####################################################################################
#
# "crawls" a single directory going thru its files looking for ready_filename
#
####################################################################################
def parse_directory(arg, dirname, names):
	try:
		# --- does the ready_to_run exist?
		cmd = "ls " + dirname + "/" + ready_filename;

		# --- execute!a
		if (os.system(cmd) == 0):
			process_pardir(dirname);
	except:
		e = sys.exc_info()[0];
		print "failed to process directory " + dirname + " message: ", e;


####################################################################################
#
# The actual directory crawler
#
####################################################################################

# --- go thru each directory at top
os.path.walk(root_search_path, parse_directory, "");


