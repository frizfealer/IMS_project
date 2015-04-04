# IMS_project
## Before using this package:
There are several depend libraries in R and Matlab that are needed. Please install these packages.
* R: 
  1. MALDIquant
  2. MALDIquantForeign
  3. MassSpecWavelet
* MATLAB:
  - For **MOLDL**
    1. [glmnet](http://web.stanford.edu/~hastie/glmnet_matlab/) 
    2. [SLEP](http://www.public.asu.edu/~jye02/Software/SLEP/overview.htm)
	3. [minFunc](http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html)
    4. [Ipopt](https://projects.coin-or.org/Ipopt/wiki/MatlabInterface)
  - For sparse-NNMF
    1. [The NMF MATLAB Toolbox] (https://sites.google.com/site/nmftool/)

## Using this package:
There are several folders in this project that are relevant to the work in ISMB2015.
* Bruker_files_conversion
* dictionary_learning_ADMM
* dictionary_learning_analysis
* dictionary_learning_syntheticData
* example
* testFunctions
* utility

# Input data format:
  The input data should be of the mzML format. This is an open XML-based format so our program can process it. 
  To covert from the Bruker proprietary file format to the mzML, One can use [msconvert](http://proteowizard.sourceforge.net/tools.shtml). The details of how to do this is in the website.
# Preprocessing:
  1. File reading and peak picking:
  after the data is converted into mzML format, using the R function MALDI_IMS_preprocessing to do the file reading and peak picking. 
  There will be four csv files in this folder.
  2. Files aggreation and binning: 
  The files in the output folder are aggregated into a file using the MATLAB function preprocess2dataCube for further MATLAB processing.
  The processed file becomes a mat file. Later on this file is loaded to bin with the MATLAB function binningDataCube.
  The example of using it can be seem in the file expScript_20150401.m under the folder "example\realExp_100831_348_136_12_40hr_0-1XLB_LN" 
  or "example\realExp_100831_348_136_12_40hr_0-1XLB_LP".
  
  