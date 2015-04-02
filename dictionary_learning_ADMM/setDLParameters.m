function [ param ] = setDLParameters()
%--------------------------------------------------------------------------
% binningDataCube: set dictionary learning parameters to default values
%--------------------------------------------------------------------------
% DESCRIPTION:

% INPUT ARGUMENTS:
% OUTPUT ARGUMENTS:
%   param data structre, has fields of:
%       OUTER_IT_NUM: # outer loop of dictionary learning
%       ADMM_IT_NUM: # loop of ADMM update on W
%       UP_D_IT_NUM: # loop of update on D
%       HES_FLAG: # wheter to use Hessian on updating D, currently deprecated
%           CLUSTER_NUM: # CPU used in dictionary learning, currently
%           deprecated, always set as one.
%           CLUSTER_NAME: name of cluster, currently deprecated.
%       SAVE_TEM_PERIOD: # of loops to save a snapshot of curret learned
%       var.
%           INNTER_IT_RED_FLAG: A flag of whether to have less loops in update
%           W and update D, currently deprecated.
%           LATE_UPDATE_FLAG: A flag of whether doing late update, currently
%           deprecated.
%           LATE_UPDATE_PERCENT: the percentage of dictionary elements used in
%           the late update situation, currently deprecated.
%           LATE_UPDATE_INTRES: the threshold values of dictionary elements
%           used in the late update situation, currently deprecated.
%       USE_L1_FLAG: a flag of whether to use L1 in fused W
%       INIT_D: initiailzation method for D
%       D_HIST_PATH: the file path to save D in each loops, if set [], it 
%       means no saving.
%       W_HIST_PATH: the file path to save W in each loops, if set [], it 
%       means no saving.
%           newWInfo: newW for the second time update. Under construction
%       W_LOWER_BOUND
%       D_LOWER_BOUND
%       LP_TOL
%       W_TOL
%       D_TOL
%           linkFunc: link function in the model. Under construction

param = [];
param.OUTER_IT_NUM = 100;
param.ADMM_IT_NUM = 100;
param.UP_D_IT_NUM = 200;
param.HES_FLAG = 0;
param.CLUSTER_NUM = 1;
param.CLUSTER_NAME = 'local';
param.SAVE_TEM_PERIOD = 1;
param.INNER_IT_RED_FLAG = 0;
param.LATE_UPDATE_FLAG = 0;
param.LATE_UPDATE_PERCENT = 0.2;
param.LATE_UPDATE_INTTHRES = 0.8;
param.USE_L1_FLAG = 0;
param.INIT_D = 'NNMF';
param.D_HIST_PATH  = [];
param.W_HIST_PATH = [];
%param.D_ION_NAME = [];
param.W_LOWER_BOUND = 1e-2;
param.D_LOWER_BOUND = 1e-2;
param.newWInfo = [];
param.LP_TOL = 1e-6;
param.W_TOL = 5*1e-3;
param.D_TOL = 5*1e-3;
param.linkFunc = 'identity';
%param.D_CONSTRAINTS = 'L1';
%param.OFFSET_FLAG = 1;
end

