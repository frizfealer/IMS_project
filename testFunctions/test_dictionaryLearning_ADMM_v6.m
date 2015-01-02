function [ output_args ] = test_dictionaryLearning_ADMM_v6( testCaseNum )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
if testCaseNum == 1
    load( 'D:\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LN\100831_348_136_12_40hr_0_1XLB_LN_input.mat' );
    load( 'D:\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LN\exp_100831_348_136_12_40hr_0_1XLB_LN_4_res.mat' );
    [ param ] = setDLParameters();
    lambda = 1e-32; theta = 1e-32; phi = 1e-32;
    param.ADMM_IT_NUM = 3;
    param.UP_D_IT_NUM = 30;
    [ expRec ] = dictionaryLearning_ADMM_v6( dataCube, [], nDTemplate, [], lambda, theta, phi, aMatrix, BlkDS, [], 'D:\temp.mat', param );
end
end

