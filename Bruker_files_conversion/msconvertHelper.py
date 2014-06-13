import fnmatch
import os
import sys
import subprocess
import re

inPath = str(sys.argv[1]);
outPath = str(sys.argv[2]);
opt = int(sys.argv[3]);
print 'Input path folder is: %s' % (inPath)
if opt == 1:
	print 'rename from 1SLin -> 1SRef'
elif opt ==2:
	print 'rename from 1SRef -> 1SLin'

print 'Get all file path needed to rename...'
matches = []
#Flag = False
if opt == 1:
	for root, dirnames, filenames in os.walk(inPath):
		#print os.path.join(root)	
		for dirname in fnmatch.filter(dirnames, '1SLin'):
			print os.path.join(root, dirname)
			oName = os.path.join(root, dirname)
			nName = os.path.join(root, '1SRef' )
			os.rename( oName, nName )
			#Flag = True
			#break
		#if Flag == True:
		#	break
elif opt == 2:
	for root, dirnames, filenames in os.walk(inPath):
		for dirname in fnmatch.filter(dirnames, '1SRef'):
			print os.path.join(root, dirname)
			oName = os.path.join(root, dirname)
			nName = os.path.join(root, '1SLin' )
			os.rename( oName, nName )
elif opt == 3:
	for root, dirnames, filenames in os.walk(inPath):
		for dirname in fnmatch.filter(dirnames, '1SLin'):
			print os.path.join(root, dirname)
			oName = os.path.join(root, dirname)
			nName = os.path.join(root, '1SRef' )
			os.rename( oName, nName )
	minPath = '"'+inPath+'"'
	outPath = '"'+outPath+'"'
	cmd = 'C:\Program Files\ProteoWizard\ProteoWizard 3.0.6212\msconvert ' + minPath + ' -o ' + outPath + ' --mzML --filter "peakPicking true 1-"' 
	print cmd
	#FNULL = open(os.devnull, 'w')
	subprocess.call(cmd, shell=False)
	for root, dirnames, filenames in os.walk(inPath):
		for dirname in fnmatch.filter(dirnames, '1SRef'):
			print os.path.join(root, dirname)
			oName = os.path.join(root, dirname)
			nName = os.path.join(root, '1SLin' )
			os.rename( oName, nName )
