function [ WAdmm_class ] = initWADMM( InputData, mLen)
%--------------------------------------------------------------------------
%initWADMM: initialize the data structure of W ADMM
%--------------------------------------------------------------------------
% DESCRIPTION:
%
% INPUT ARGUMENTS: 
%   InputData, the InputData data structure.
%   mLen, # molecules in the dictionary.
% OUTPUT ARGUMENTS: 
%   WAdmm_class, a data structure has parameters listed below.
%       itNum: # loop of ADMM update on W
%       lambda: a hyperparameter for controling W's sparsity
%       theta: a hyperparameter for controling W's fused terms.
%       W_TOL: a terminating criterion for updating. if the difference
%       between W before updating and W after updating is smaller then
%       maxW*W_TOL, then break.
%       D_LOWER_BOUND: the lower bound of a dictionary element to be
%       considered as a valid element. If max(D_element) < D_LOWER_BOUND,
%       it does not consider this element.
%       LINK_FUNC: link function in the model. Default is 'identity'. 
%       Others are under construction.
%       EPS_ABS_PRI, EPS_REL_PRI: epsilon constants for primal variables.
%       EPS_ABS_DUAL, EPS_REL_DUAL: epsilon constants for dual variables.
%       L1_FLAG: a flag of whether to use L1 in fused W
%       RHO_INIT: the itial rho for the update.
%       stuckFlag: a flag indicate whether the process stuck while
%       updating. "Stuck" means cannot make LP smaller within the
%       2*WAdmm_class.itNum iterations.
%   WAdmm_class also has model data fields listed below.
%       z0, z1, z2, W, W0

%% pour out input data
Y = InputData.dataCube;
BlkDS = InputData.BlkDS;
[~, hei, wid] = size( Y );
aMatrix = InputData.aMatrix;
if ~isempty(InputData.Rblk)
    Rblk = InputData.Rblk;
else %for test data
    Rblk = [];
end

%% initialize parameters
WAdmm_class.itNum= 100;
WAdmm_class.lambda =0;
WAdmm_class.theta = 0;
WAdmm_class.WHistFlag = 0;
WAdmm_class.W_TOL = 1e-3;
WAdmm_class.D_LOWER_BOUND = 1e-2;
WAdmm_class.LINK_FUNC = 'identity';
WAdmm_class.EPS_ABS_PRI = 1e-4;
WAdmm_class.EPS_REL_PRI = 1e-4;
WAdmm_class.EPS_ABS_DUAL = 1e-4;
WAdmm_class.EPS_REL_DUAL = 1e-2;
WAdmm_class.L1Flag = 0;
WAdmm_class.RHO_INIT = 1e-6;
WAdmm_class.stuckFlag = 0;
%% initialize models
%setting W, W0, z0, z1, z2
if strcmp( WAdmm_class.LINK_FUNC, 'log' ) == 1
    WAdmm_class.z0 = log(Y);
elseif strcmp( WAdmm_class.LINK_FUNC, 'identity' ) == 1
    WAdmm_class.z0 = Y;
    %if z0 is in the colony's area (i) but = 0, set the m/z channels in
    %z0(i) to 1
    for i = 1:length(BlkDS.indMap(:))
        if BlkDS.indMap(i) == 1
            idx =  WAdmm_class.z0(:, i) == 0 ;
            WAdmm_class.z0(idx, i) = 1;
        end
    end
elseif strcmp( WAdmm_class.LINK_FUNC, 'log_gaussain' ) == 1
    WAdmm_class.z0 = log(Y);
elseif strcmp( WAdmm_class.LINK_FUNC, 'negative_binomial' ) == 1
    WAdmm_class.z0 = log(Y);
end
%if not in traing set, we shoud not initialize z0 according to it.
%take the average values
idx = intersect( find(aMatrix==1), find(BlkDS.indMap==1) );
tmp = sum( WAdmm_class.z0(:, idx), 2 ) / ( length( idx ) );
for i = 1:hei*wid
    %if is in data area, but not in training set,
    %set it to average spectrum
    if BlkDS.indMap(i) == 1 && aMatrix(i) == 0
        WAdmm_class.z0(:, i) = tmp;
    end
end
WAdmm_class.z1 = zeros( mLen, wid*hei );
if ~isempty( Rblk )
    WAdmm_class.z2 = cell( BlkDS.blkNum, 1 );
    for i = 1:BlkDS.blkNum
        WAdmm_class.z2{i} = zeros( length(Rblk{i}), mLen );
    end
else
    WAdmm_class.z2 = [];
end
WAdmm_class.W = zeros( mLen, wid*hei );
WAdmm_class.W0 = ones( hei, wid );
%setting WHistCell
WAdmm_class.WHistCell = cell( WAdmm_class.itNum, 1 );

end

