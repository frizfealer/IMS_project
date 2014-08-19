import fnmatch
import os
import sys
import subprocess
import re

inPath = str(sys.argv[1]);
outPath = str(sys.argv[2]);
print 'Input path folder is: %s' % (inPath)
print 'Output Path folder is: %s' % (outPath)

if not os.path.exists(outPath):
    os.makedirs(outPath)

	
print 'Get all fid file path...'
filePaths = os.listdir( inPath );

#for R script format
outPath = outPath.replace( '\\', '\\\\' );
for filePath in filePaths:
	realPath = os.path.join(inPath, filePath);
	#for R script format
	realPath = realPath.replace( '\\', '\\\\' );
	cmd = 'D:\\Users\\yharn\\Documents\\R\\R-3.1.1\\bin\\x64\\Rscript D:\\usePreprocessing.R "' + realPath + '" "' + outPath + '\\\\"';
	print cmd
	#FNULL = open(os.devnull, 'w')
	subprocess.call(cmd, shell=False)


# f = open('fid_path.txt', 'w')
# for item in matches:
	# f.write("%s\n" % item)
# f.close()

# for item in matches:
	# print "processing %s" % (item)
	# m = re.search( r'.*\\(.*)\\1\\1SLin\\fid', item );
	# outfileName = outPath +'\\'+m.group(1)+'.mzML'
	# f = open( outfileName, 'w')
	# f.close()
	# inFile = '"' + item + '"'
	# outPathC = '"' + outPath + '"'
	# cmd = 'C:\Program Files\ProteoWizard\ProteoWizard 3.0.6212\msconvert ' + inFile + ' -o ' + outPathC + ' --outfile ' + m.group(1)+'.mzML'
	# #print cmd
	# FNULL = open(os.devnull, 'w')
	# subprocess.call(cmd, stdout=FNULL, stderr=FNULL, shell=False)
# #FNULL = open(os.devnull, 'w')    #use this if you want to suppress output to stdout from the subprocess
# filename = "my_file.dat"
# args = "RegressionSystem.exe -config " + filename
# subprocess.call(cmd, shell=False)
# subprocess.call(cmd, stdout=FNULL, stderr=FNULL, shell=False)