function [ output_args ] = realExperiment_100827_WP13_WP17_2( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
load( '../100831_113_333_136_26hr_0_1XLB_2_inputs_env.mat' );
addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/SLEP_package_4.1/SLEP' ) );
addpath( '~/Ipopt_3118' );
param.CLUSTER_NAME = 'killdevil1024'; %usually this is the default name
param.CLUSTER_NUM = 64; %number of cluster uses, usually 12
snapPath = 'realExperiment_100827_WP13_WP17_2_temp.mat';
param.OUTER_IT_NUM = 100; %number of  outer loop
elambda =  3.7486; etheta =  4.8046; ephi =  4.1983;
[ expRecReal ] = dictionaryLearning_ADMM_v4( tDataCube, [], tDTemplate, [], elambda, etheta, ephi, aMatrix, tBlkDS, [], snapPath, param );
save( 'realExperiment_100827_WP13_WP17_2_res.mat', 'expRecReal', '-v7.3' ) ;

end

