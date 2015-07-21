# IMS_project
## Before using this package:
There are several depend libraries in R and Matlab that are needed. Please install these packages.
* R: 
  1. MALDIquant
  2. MALDIquantForeign
  3. MassSpecWavelet
  4. mzR
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

## Preprocessing:
  1. The input data, whaterever intruments of IMS is using, should be converted into the mzML format. This is an open XML-based format so our program can process it. 
  To covert from the Bruker proprietary file format to the mzML, One can use [msconvert](http://proteowizard.sourceforge.net/tools.shtml). The details of how to do this is in Bruker_files_conversion\commands.txt.
  If you need the mzML file in the example, you can download them from: https://drive.google.com/open?id=0B06g-wOYKgZ8ald2VGtKUnJIYVE and https://drive.google.com/open?id=0B06g-wOYKgZ8amtnRmNLRnE1MDA.
  2. Then we use Bruker_files_conversion\usePreprocessing.R to read mzML file and do peak picking. Again, the details of how to do this is in Bruker_files_conversion\commands.txt.
  3. The steps afterwards are shown in exp_100831_348_136_12_40hr_0_1XLB_LN_basicSettings_v2.m and exp_100831_348_136_12_40hr_0_1XLB_LP_basicSettings_v2.m

## Main program:
  The example of using it can be seem in the file **expScript_20150401.m** under the folder "example\realExp_100831_348_136_12_40hr_0-1XLB_LN" 
  or "example\realExp_100831_348_136_12_40hr_0-1XLB_LP".
## 
  
  