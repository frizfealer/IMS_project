function [ output_args ] = realExperiment_2013_Bmyc_Paeni_Early_LP_3(  input_args )

%% run on different hyper-parameters limits
load( '../2013_Bmyc_Paeni_Early_LP_2_inputs_env.mat' );
addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/SLEP_package_4.1/SLEP' ) );
param.CLUSTER_NAME = 'killdevil1024'; %usually this is the default name
param.CLUSTER_NUM = 64; %number of cluster uses, usually 12
%% generate BlkDS, a data structure for bacteria community location information
%BlkDS has the 4 fields
%blkNum: the bacteria community number 
%B2GMap: a mapping from bacteria community to grid
%G2BMap: a mapping from grid to bacteria community
%indMap: a logical matrix, as the same size of the grid, with 1 means the
%grids having signals, 0 means empty
[SLEN, IHEIGHT, IWIDTH] = size( dataCube );
[ BlkDS ] = conBLKDS( dataCube );
elambda = 1.9303; etheta = 2.0707; ephi = 4.2295;
snapPath = 'realExperiment_2013_Bmyc_Paeni_Early_LP_3_temp.mat';
aMatrix = ones( IHEIGHT, IWIDTH );
param.OUTER_IT_NUM = 100; %number of  outer loop
[ expRecReal ] = dictionaryLearning_ADMM_v4( dataCube, [], pDTemplate, [], elambda, etheta, ephi, aMatrix, BlkDS, [], snapPath, param );
save( 'realExperiment_2013_Bmyc_Paeni_Early_LP_3_res.mat', 'expRecReal', '-v7.3' ) ;