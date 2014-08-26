function [ output_args ] = realExperiment_100827_WP13_WP17_2( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
load( 'realExperiment_100827_WP13_WP17_1_results_env.mat' );
addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/SLEP_package_4.1/SLEP' ) );
addpath( '~/Ipopt_3118' );
param.CLUSTER_NAME = 'killdevil32'; %usually this is the default name
param.CLUSTER_NUM = 32; %number of cluster uses, usually 12
snapPath = 'realExperiment_100827_WP13_WP17_3_temp.mat';
param.OUTER_IT_NUM = 100; %number of  outer loop
param.UP_D_IT_NUM = 200;
%for killdevil setup
[ expRecReal ] = dictionaryLearning_ADMM_v4( tDataCube, [], tDTemplate, [], elambda, etheta, ephi, aMatrix, tBlkDS, [], snapPath, param );
save( 'realExperiment_100827_WP13_WP17_3_res.mat', 'expRecReal', '-v7.3' ) ;

end

