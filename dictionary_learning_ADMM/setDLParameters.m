function [ param ] = setDLParameters()
%setDLParameters set dictionary learning parameters to default values
param = [];
param.OUTER_IT_NUM = 100;
param.ADMM_IT_NUM = 200;
param.UP_D_IT_NUM = 200;
param.HES_FLAG = 1;
param.CLUSTER_NUM = 1;
param.CLUSTER_NAME = 'local';
param.SAVE_TEM_PERIOD = 1;
param.INNER_IT_RED_FLAG = 0;
param.LATE_UPDATE_FLAG = 0;
param.LATE_UPDATE_PERCENT = 0.2;
param.LATE_UPDATE_INTTHRES = 0.8;
param.USE_L1_FLAG = 1;
param.INIT_D = 'NNMF';
param.D_HIST_PATH  = [];
param.W_HIST_PATH = [];
param.D_ION_NAME = [];
param.W_LOWER_BOUND = 1e-2;
param.D_LOWER_BOUND = 1e-2;
param.newWInfo = [];
param.LP_TOL = 1e-6;
param.W_TOL = 1e-3;
param.D_TOL = 1e-3;
end

