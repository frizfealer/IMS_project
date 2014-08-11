function [ output_args ] = realExperiment_100831_113_333_136_26hr_0_1XLB_3( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
load( '../100831_113_333_136_26hr_0_1XLB_2_inputs_env.mat' );
addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/SLEP_package_4.1/SLEP' ) );
param.CLUSTER_NAME = 'killdevil1024'; %usually this is the default name
param.CLUSTER_NUM = 64; %number of cluster uses, usually 12
snapPath = 'realExperiment_100831_113_333_136_26hr_0_1XLB_3_temp.mat';
[SLEN, IHEIGHT, IWIDTH] = size( tDataCube );
[ BlkDS ] = conBLKDS( tDataCube );
aMatrix = ones( size( tDataCube, 2 ), size( tDataCube, 3) );
param.OUTER_IT_NUM = 100; %number of  outer loop
elambda =  3.7486; etheta =  4.8046; ephi =  4.1983;
[ expRecReal ] = dictionaryLearning_ADMM_v4( tDataCube, [], tDTemplate, [], elambda, etheta, ephi, aMatrix, BlkDS, [], snapPath, param );
save( 'realExperiment_100831_113_333_136_26hr_0_1XLB_3_res.mat', 'expRecReal', '-v7.3' ) ;

end

