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
  The input data should be of the mzML format. This is an open XML-based format so our program can process it. To covert from a proprietary file format to the mzML, One can use [msconvert](http://proteowizard.sourceforge.net/tools.shtml). The details of how to do this is in the website.
# Preprocessing:
  1. File reading and peak picking:
    After the data is converted into mzML format, using the R function MALDI_IMS_preprocessing( iPath, outputPath, 0, 0 ) to do the file reading and peak picking. The iPath is the file path of the input data from the step 1. Rhe outputPath shoud be the path of output folder (e.g. D:\\result\\). There will be four csv files in this folder.