function [ expRec ] = dictionaryLearning_ADMM_v6( inY, initVar, DTemplate, logFY, lambda, theta, phi, aMatrix, BlkDS, filePath, snapFilePath, param, varargin )
%dictionaryLearning_ADMM_v6 dictionary learning with ADMM formulation
%with set z0 = D*W+W0, z1 = w1, z2 = RW';
%Dictioanry Learning with Poisson distribution
%input:
% inY [s h w]
% initVar data structure with all initialization, has the fields:
% initVar.outW initVar.outW0 initVar.outD 
% initVar.u0 initVar.u1 initVar.z0 initVar.z1, initVar.z2, initVar.u2
%DTemplate the dictionary pattern you want to have [s m]
%We choose the appropriate dictionary pattern before this function
%lambda, theta, phi: scalar
%filePath: the output file path, if output to screen, set []
%aMatrix [h w] a Matrix indicate which grid is held out, 1 means trainning,
%0 means testing
%snapFilePath: file path for snapshot according to SAVE_TEM_PERIOD
%BlkDS, structures with two fields,
%   bIdx is the block to location index
%   mapW2B is the location to block index
%   blkNum the number of the blocks
%output:
%expRec is a structure contains 
%LPAry [itNum+1, 1], a log posterior record after each iteration, the first
%one is the log posterior from orignal iteration
%outW, outW0, outD, z0, z1, u0, u1
%rho, theta, lambda
[~, hei, wid] = size(inY);
[mLen] = size(DTemplate ,2);

%% set output path, and iteration number
if ~isempty(filePath)
    fID = fopen(filePath, 'w');
else
    fID = 1;
end

if isfield( param, 'OUTER_IT_NUM')
    OUTER_IT_NUM = param.OUTER_IT_NUM;
else
    OUTER_IT_NUM = 100;
end
if isfield( param, 'ADMM_IT_NUM' )
    ADMM_IT_NUM = param.ADMM_IT_NUM;
else
    ADMM_IT_NUM = 200;
end
if isfield( param, 'UP_D_IT_NUM' )
    UP_D_IT_NUM = param.UP_D_IT_NUM;
else
    UP_D_IT_NUM = 200;
end
if isfield( param, 'HES_FLAG' )
    HES_FLAG = param.HES_FLAG;
else
    HES_FLAG = 1;
end
if isfield( param, 'CLUSTER_NUM' )
    CLUSTER_NUM = param.CLUSTER_NUM;
else
    CLUSTER_NUM = 1;
end
if isfield( param, 'SAVE_TEM_PERIOD' )
    SAVE_TEM_PERIOD = param.SAVE_TEM_PERIOD;
else
    SAVE_TEM_PERIOD = 3;
end
if isfield( param, 'INNER_IT_RED_FLAG' )
    INNER_IT_RED_FLAG = param.INNER_IT_RED_FLAG;
else
    INNER_IT_RED_FLAG = 0;
end
% if isfield( param, 'LATE_UPDATE_FLAG' )
%     LATE_UPDATE_FLAG = param.LATE_UPDATE_FLAG;
% else
%     LATE_UPDATE_FLAG = 0;
% end
% if isfield( param, 'LATE_UPDATE_PERCENT' )
%     LATE_UPDATE_PERCENT = param.LATE_UPDATE_PERCENT;
% else
%     LATE_UPDATE_PERCENT = 0.2;
% end
% if isfield( param, 'LATE_UPDATE_INTTHRES' )
%     LATE_UPDATE_INTTHRES = param.LATE_UPDATE_INTTHRES;
% else
%     LATE_UPDATE_INTTHRES = 0.8;
% end
if isfield( param, 'CLUSTER_NAME' )
    CLUSTER_NAME = param.CLUSTER_NAME;
else
    CLUSTER_NAME = 'local';
end
if isfield( param, 'INIT_D' )
    INIT_METHOD = param.INIT_D;
else
    INIT_METHOD = 'NNMF';
end
if isfield( param, 'D_HIST_PATH' )
    D_HIST_PATH = param.D_HIST_PATH;
else
    D_HIST_PATH = [];
end
if isfield( param, 'W_HIST_PATH' )
    W_HIST_PATH = param.W_HIST_PATH;
else
    W_HIST_PATH = [];
end
DHistCell = cell( OUTER_IT_NUM+1, 1 );
WHistCell = cell( OUTER_IT_NUM+1, 1 );
if isfield( param, 'D_ION_NAME' )
    D_ION_NAME = param.D_ION_NAME;
else
    D_ION_NAME = [];
end
if isfield( param, 'LATE_UPDATE_IT_PREC' )
    LATE_UPDATE_IT_PREC = param.LATE_UPDATE_IT_PREC;
else
    LATE_UPDATE_IT_PREC = 0.8;
end
if isfield( param, 'USE_L1_FLAG' )
    USE_L1_FLAG = param.USE_L1_FLAG;
else
    USE_L1_FLAG = 1;
end
if isfield( param, 'W_lOWERBOUND' )
    W_LOWER_BOUND = param.W_LOWER_BOUND;
else
    W_LOWER_BOUND = 1e-2;
end
if isfield( param, 'LP_TOL' )
    LP_TOL = param.LP_TOL;
else
    LP_TOL = 1e-6;
end
if isfield( param, 'W_TOL' )
    W_TOL = param.W_TOL;
else
    W_TOL = 1e-3;
end
if isfield( param, 'D_TOL' )
    D_TOL = param.D_TOL;
else
    D_TOL = 1e-3;
end

LPAry = zeros( 1, OUTER_IT_NUM+1 );

%rhoCell will consist of MAXIT iteration number of rhoAry
% rhoCell = cell( OUTER_IT_NUM, 1 );
%residual Record Cell will consist of MAXIT iteration number of resRecAry
% resRecCell = cell(OUTER_IT_NUM, 1);

%% initialize all variable
if isempty(initVar)
    [ D ] = initD( inY, DTemplate, INIT_METHOD, D_ION_NAME );
    DHistCell{1} = D;
    W = sparse( zeros( mLen, hei*wid ) );
    W0 = sparse( zeros( hei, wid ) );
    WHistCell{1} = sparse( zeros( mLen, hei*wid ) );
else
    D = initVar.D;
    DHistCell{1} = D;
    W = initVar.W;
    WHistCell{1} = W;
    W0 = initVar.W0;
    z0 = initVar.z0; z1 = initVar.z1; z2 = initVar.z2;
end

%% managing matlab pool
if matlabpool('size') == 0 && CLUSTER_NUM > 1
    matlabpool( 'open', CLUSTER_NAME, CLUSTER_NUM );
end
%%managing late update for some dictionary elements
% if LATE_UPDATE_FLAG == 1
%     [ rSet ] = estimateDRelevant( inY, outD, DTemplate, BlkDS.indMap, LATE_UPDATE_PERCENT, LATE_UPDATE_INTTHRES, D_ION_NAME  );
% end
%% main program start
MEAN_FLAG = 0;
%compute scaleFactor
tmp = log(inY); tmp(tmp==-inf)=0; 
idx = intersect( find(aMatrix==1), find(BlkDS.indMap==1) );
tmp2 = sum( tmp(:, idx), 2 ) / ( length( idx ) );
for i = 1:hei*wid
    if BlkDS.indMap(i) == 1 && aMatrix(i) == 0
        tmp(:, i) = tmp2;
    end
end
scaleFactor =  1 / ( max( inY(:) ) / max( tmp(:) ) );
LPAry(1) = LP_DL_Poiss( aMatrix, inY, W, W0, D, lambda, phi, theta, scaleFactor, logFY, MEAN_FLAG );
fprintf( 'parameters: outer iteration number = %d, ', OUTER_IT_NUM );
fprintf( 'ADMM iteration number = %d, D-update iteration number = %d, Hessian flag for update D = %d, ', ADMM_IT_NUM, UP_D_IT_NUM, HES_FLAG );
fprintf( 'Cluster number used to update = %d, ', CLUSTER_NUM );
fprintf( 'SAVE_TEM_PERIOUD = %d, Inner iteration reduce flag = %d, ', SAVE_TEM_PERIOD, INNER_IT_RED_FLAG );
%fprintf( 'LATE_UPDATE_FLAG = %d, LATE_UPDATE_PERCENT=%g, LATE_UPDATE_INTTHRES = %g, ', LATE_UPDATE_FLAG, LATE_UPDATE_PERCENT, LATE_UPDATE_INTTHRES  );
fprintf( ['INIT_D= ', INIT_METHOD, '.\n'] );
for it = 1:OUTER_IT_NUM
    fprintf( fID, 'iteration %d\n', it );
    if INNER_IT_RED_FLAG == 1
        if it < floor(OUTER_IT_NUM*LATE_UPDATE_IT_PREC)
            M_ADMM_IT_NUM = round( ADMM_IT_NUM / 2 );
            M_UP_D_IT_NUM = round( UP_D_IT_NUM / 2 );
        else
            M_ADMM_IT_NUM = ADMM_IT_NUM;
            M_UP_D_IT_NUM = UP_D_IT_NUM;
        end
    else
        M_ADMM_IT_NUM = ADMM_IT_NUM;
        M_UP_D_IT_NUM = UP_D_IT_NUM;
    end
    prevW = W; prevW0 = W0; prevD = D;
    %% update W
    if it == 1 
        if ~isempty(initVar)
            curVar.W = initVar.W; curVar.W0 = initVar.W0;
            curVar.z0 = initVar.z0; curVar.z1 = initVar.z1; curVar.z2 = initVar.z2;
        else
            curVar = [];
        end
    else
        curVar.W = uW_Res.W; curVar.W0 = uW_Res.W0;
        curVar.z0 = uW_Res.z0; curVar.z1 = uW_Res.z1; curVar.z2 = uW_Res.z2;
    end
    %the second to the last parameter is a flag for each W output in ADMM steps
    %the last parameter is the tolerance of w in ADMM steps
    uW_Res = updateW_ADMM_v3( inY, D, aMatrix, M_ADMM_IT_NUM, lambda, theta, USE_L1_FLAG, logFY, curVar, scaleFactor, 0, 1e-2 );
    W = uW_Res.W;
    W0 = uW_Res.W0;
    z0 = uW_Res.z0; z1 = uW_Res.z1; z2 = uW_Res.z2;
    
    %% if W is too low, make corresponding D to zero
    for i = 1:size(W, 1)
        if max(W(i,:)) < W_LOWER_BOUND
            D( DTemplate(:,i)==1, i) = 1e-32;
        end
    end
    %% update D
    validMap = BlkDS.indMap .* aMatrix;
    [ D ]= updateD_v8_ipopt( inY, W, W0, D, DTemplate, validMap, HES_FLAG, phi, scaleFactor, M_UP_D_IT_NUM, W_LOWER_BOUND );

    LPAry(it+1) = LP_DL_Poiss( aMatrix, inY, W, W0, D, lambda, phi, theta, scaleFactor, logFY, MEAN_FLAG );
    tmp1 = max( abs( W(:)-prevW(:) ) );
    tmp2 = max( abs( W0(:)-prevW0(:) ) );
    tmp3 = max( abs( D(:)-prevD(:) ) );
    if mod(it, SAVE_TEM_PERIOD) == 0
        expRec.D = D; expRec.W = W; expRec.W0 = W0;
        expRec.LPAry = LPAry; expRec.z0 = z0; expRec.z1 = z1; expRec.z2 = z2;
        expRec.theta = theta; expRec.lambda = lambda; expRec.phi = phi;
        expRec.diffW = tmp1; expRec.diffW0 = tmp2; expRec.diffD = tmp3;
        expRec.param = param;
%         expRec.rhoCell = rhoCell; expRec.resRecCell = resRecCell;
        save( snapFilePath, 'expRec' );
    end
    if ~isempty( D_HIST_PATH )
        DHistCell{it+1} = D;
        save( D_HIST_PATH, 'DHistCell' );
    end
    if ~isempty( W_HIST_PATH )
        WHistCell{it+1} = W;
        save( W_HIST_PATH, 'WHistCell', '-v7.3' );
    end
    % ALL conditions must be satisfied to stop outer loops
    if abs(LPAry(it+1)-LPAry(it)) <= LP_TOL && ...
        ( max( tmp1, tmp2 ) < W_TOL &&  tmp3 < D_TOL )
        break;
    end
end
if ~isempty(filePath)
    fclose(fID);
end
expRec.D = D; expRec.W = W; expRec.W0 = W0;
expRec.LPAry = LPAry; expRec.z0 = z0; expRec.z1 = z1; expRec.z2 = z2;
expRec.theta = theta; expRec.lambda = lambda; expRec.phi = phi; 
expRec.diffW = tmp1; expRec.diffW0 = tmp2; expRec.diffD = tmp3;
expRec.param = param;
%expRec.rhoCell = rhoCell; expRec.resRecCell = resRecCell;

end

