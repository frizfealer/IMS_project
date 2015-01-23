import fnmatch
import os
import sys
import subprocess
import re
import time
#inPath = str(sys.argv[1]);
#print 'Input path folder is: %s' % (inPath)
#os.chdir(inPath)
for root, dirnames, filenames in os.walk(os.getcwd()):
	print 'Get all files in this folder...'
[tmp, folderName] = os.path.split(os.getcwd())
print folderName
logFolderName = folderName+'_logs'
logFNPath = '../'+logFolderName+'/'
if not os.path.exists(logFNPath):
    os.makedirs(logFNPath)
for file in filenames:
	#cmd = 'matlab -nodesktop -nodisplay -nosplash -r "' + file[:-2] + '; exit" -logfile ' + file[:-2]+'.log &'
	cmd = 'bsub matlab -q day -nodisplay -nosplash -singleCompThread -r ' + file[:-2] + ' -logfile ' + logFNPath + file[:-2]+'.log';

	print cmd
	#subprocess.call(cmd, shell=False)
	#subprocess.Popen(cmd, cwd = inPath);
	#os.system( cmd )
	time.sleep(1)
