function [ expRec ] = dictionaryLearning_ADMM_v6_2( inY, initVar, DTemplate, logFY, lambda, theta, phi, aMatrix, BlkDS, filePath, snapFilePath, param )
%--------------------------------------------------------------------------
% dictionaryLearning_ADMM_v6_2: dictionary learning with ADMM formulation 
% version 6 for ISMB 2015
%--------------------------------------------------------------------------
% DESCRIPTION:
%   The main function of the algorithm. Poisson distribution.
%   ADMM setting: z0 = D*W+W0, z1 = w1, z2 = RW';
% INPUT ARGUMENTS:
%   inY, input data cube, size of [s w h], that is #(m/z)*width*height
%   initVar, data structure with all initialization, has the fields:
%   initVar.outW initVar.outW0 initVar.outD 
%   initVar.u0 initVar.u1 initVar.z0 initVar.z1, initVar.z2, initVar.u2
%   DTemplate, the dictionary pattern you want to have [s m], generatd from
%   genDTemplate_v3
%   lambda, theta, phi: scalar, hyper-parameters
%   filePath, the output file path, if output to screen, set []
%   aMatrix, [h w] a indicator matrix generated from leaveoutCBPatternsData_v2
%   snapFilePath, file path for snapshot according to SAVE_TEM_PERIOD
%   BlkDS, data cube structre variable generated from conBLKDS
%   param, all parameters
% OUTPUT ARGUMENTS:
%   expRec, a data structure contains:
%   1. D, dictionary. 2. W, weight. 3. W0, intercept.
%   4. z0, z1, z2, please see the formulation of ADMM
%   5. lambda, theta, phi, hyper-parameters used in this model.
%   6. diffW: maximun difference of W between last two loops.
%   7. diffD: maximun differece of D between last two loops.
%   8. param: input parameters. 9. aMatrix: input a matrix
%   LPAry [itNum+1, 1], a log posterior record after each iteration, 
%   the first one is the log posterior from orignal iteration.

[~, hei, wid] = size(inY);
[mLen] = size(DTemplate ,2);
%% set output path, and iteration number
if ~isempty(filePath)
    fID = fopen(filePath, 'w');
else
    fID = 1;
end

%% set all parameters
if isfield( param, 'OUTER_IT_NUM')
    OUTER_IT_NUM = param.OUTER_IT_NUM;
else
    OUTER_IT_NUM = 100;
end
if isfield( param, 'ADMM_IT_NUM' )
    ADMM_IT_NUM = param.ADMM_IT_NUM;
else
    ADMM_IT_NUM = 100;
end
if isfield( param, 'UP_D_IT_NUM' )
    UP_D_IT_NUM = param.UP_D_IT_NUM;
else
    UP_D_IT_NUM = 200;
end
if isfield( param, 'HES_FLAG' ) %under construction
    HES_FLAG = param.HES_FLAG;
else
    HES_FLAG = 1;
end
% if isfield( param, 'CLUSTER_NUM' ) %deprecated
%     CLUSTER_NUM = param.CLUSTER_NUM;
% else
%     CLUSTER_NUM = 1;
% end
if isfield( param, 'SAVE_TEM_PERIOD' )
    SAVE_TEM_PERIOD = param.SAVE_TEM_PERIOD;
else
    SAVE_TEM_PERIOD = 3;
end
if isfield( param, 'INNER_IT_RED_FLAG' ) %deprecated
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
% if isfield( param, 'CLUSTER_NAME' ) %deprecated
%     CLUSTER_NAME = param.CLUSTER_NAME;
% else
%     CLUSTER_NAME = 'local';
% end
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
if isfield( param, 'D_ION_NAME' ) %deprecated
    D_ION_NAME = param.D_ION_NAME;
else
    D_ION_NAME = [];
end
if isfield( param, 'LATE_UPDATE_IT_PREC' ) %deprecated
    LATE_UPDATE_IT_PREC = param.LATE_UPDATE_IT_PREC;
else
    LATE_UPDATE_IT_PREC = 0.8;
end
if isfield( param, 'USE_L1_FLAG' ) %deprecated, always using 0
    USE_L1_FLAG = param.USE_L1_FLAG;
else
    USE_L1_FLAG = 1;
end
if isfield( param, 'W_LOWER_BOUND' )
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
if isfield( param, 'D_LOWER_BOUND' )
    D_LOWER_BOUND = param.D_LOWER_BOUND;
else
    D_LOWER_BOUND = 1e-2;
end
%for second-round dictionary learning, under construction
if isfield( param, 'newWInfo' )
    newWInfo = param.newWInfo;
else
    newWInfo = [];
end
%for link function, under construction, only valid for 'identity'
if isfield( param, 'linkFunc' )
    LINK_FUNC = param.linkFunc;
else
    LINK_FUNC = 'identity';
end
%for Dictionary element's constraints not using this for now
% if isfield( param, 'D_CONSTRAINTS' )
%     D_CONSTRAINTS = param.D_CONSTRAINTS;
% else
%     D_CONSTRAINTS = 'L1';
% end
% if isfield( param, 'OFFSET_FLAG' ) 
%     OFFSET_FLAG = param.OFFSET_FLAG;
% else
%     OFFSET_FLAG = 0;
% end


LPAry = zeros( 1, OUTER_IT_NUM+1 );
dfWVec = zeros( 1, OUTER_IT_NUM );

%rhoCell will consist of MAXIT iteration number of rhoAry
% rhoCell = cell( OUTER_IT_NUM, 1 );
%residual Record Cell will consist of MAXIT iteration number of resRecAry
% resRecCell = cell(OUTER_IT_NUM, 1);

%% initialize all variable
if isempty(initVar)
    %     if strcmp( D_CONSTRAINTS, 'L1' ) == 1
    %         [ D ] = initD( inY, DTemplate, INIT_METHOD, D_ION_NAME, 1 );
    %     else
    [ D ] = initD( inY, DTemplate, INIT_METHOD, [], 0 );
    %end
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
% if matlabpool('size') == 0 && CLUSTER_NUM > 1
%     matlabpool( 'open', CLUSTER_NAME, CLUSTER_NUM );
% end
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
% scaleFactor =  1 / ( max( inY(:) ) / max( tmp(:) ) );
if strcmp( LINK_FUNC, 'identity' ) == 1
    scaleFactor = 1;
end
% set all scaleFactor to 1
LPAry(1) = LP_DL_Poiss( LINK_FUNC, aMatrix, inY, W, W0, D, lambda, phi, theta, scaleFactor, logFY, MEAN_FLAG );
fprintf( 'parameters: outer iteration number = %d, ', OUTER_IT_NUM );
fprintf( 'ADMM iteration number = %d, D-update iteration number = %d, Hessian flag for update D = %d, ', ADMM_IT_NUM, UP_D_IT_NUM, HES_FLAG );
% fprintf( 'Cluster number used to update = %d, ', CLUSTER_NUM );
fprintf( 'SAVE_TEM_PERIOUD = %d, Inner iteration reduce flag = %d, ', SAVE_TEM_PERIOD, INNER_IT_RED_FLAG );
%fprintf( 'LATE_UPDATE_FLAG = %d, LATE_UPDATE_PERCENT=%g, LATE_UPDATE_INTTHRES = %g, ', LATE_UPDATE_FLAG, LATE_UPDATE_PERCENT, LATE_UPDATE_INTTHRES  );
fprintf( ['INIT_D= ', INIT_METHOD, '.\n'] );
fprintf( ['LINK_FUNC= ', LINK_FUNC, '.\n'] );
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
    uW_Res = updateW_ADMM_v3_2( LINK_FUNC, inY, D, aMatrix, M_ADMM_IT_NUM, lambda, theta, USE_L1_FLAG, logFY, curVar, scaleFactor, 0, W_TOL, D_LOWER_BOUND, newWInfo );
    if uW_Res.stuckFlag == 1
        expRec.WStuckFlag = 1;
    else
        expRec.WStuckFlag = 0;
    end
    %ensure sparsity in W
    uW_Res.W( uW_Res.W < W_LOWER_BOUND ) = 0;
    W = uW_Res.W;
    W0 = uW_Res.W0;
    z0 = uW_Res.z0; z1 = uW_Res.z1; z2 = uW_Res.z2;
    
    %% update D
    validMap = BlkDS.indMap .* aMatrix;
    before = LP_DL_Poiss( LINK_FUNC, aMatrix, inY, W, W0, D, lambda, phi, theta, scaleFactor, logFY, MEAN_FLAG );
    [ D ]= updateD_v8_ipopt( LINK_FUNC, 'L2_SQUARE', inY, W, W0, D, DTemplate, validMap, HES_FLAG, phi, scaleFactor, M_UP_D_IT_NUM, W_LOWER_BOUND );
    after = LP_DL_Poiss( LINK_FUNC, aMatrix, inY, W, W0, D, lambda, phi, theta, scaleFactor, logFY, MEAN_FLAG );
    
    cnt = 0;
    expRec.DStuckFlag = 0;
    while after-before > 0
        if cnt == 5
            expRec.DStuckFlag = 1;
            break;
        end
        cnt = cnt + 1;
        [ D ]= updateD_v8_ipopt( LINK_FUNC, 'L2_SQUARE', inY, W, W0, D, DTemplate, validMap, HES_FLAG, phi, scaleFactor, M_UP_D_IT_NUM, W_LOWER_BOUND );
        after = LP_DL_Poiss( LINK_FUNC, aMatrix, inY, W, W0, D, lambda, phi, theta, scaleFactor, logFY, MEAN_FLAG );
    end
    %ensure sparsity in D
    D( D < D_LOWER_BOUND ) = 0;
    
    LPAry(it+1) = LP_DL_Poiss( LINK_FUNC, aMatrix, inY, W, W0, D, lambda, phi, theta, scaleFactor, logFY, MEAN_FLAG );
    tmp1 = max( abs( W(:)-prevW(:) ) );
    tmp2 = max( abs( W0(:)-prevW0(:) ) );
    tmp3 = max( abs( D(:)-prevD(:) ) );
    if mod(it, SAVE_TEM_PERIOD) == 0
        expRec.D = D; expRec.W = W; expRec.W0 = W0;
        expRec.LPAry = LPAry; expRec.z0 = z0; expRec.z1 = z1; expRec.z2 = z2;
        expRec.theta = theta; expRec.lambda = lambda; expRec.phi = phi;
        expRec.diffW = tmp1; expRec.diffW0 = tmp2; expRec.diffD = tmp3;
        expRec.param = param; expRec.aMatrix = aMatrix;
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
    % compute df for W
%     ins = zeros( size(W, 1), 1 );
%     for i = 1:length(ins)
%         ins(i) = max( W(i,:) );
%     end
%     dfWVec(it) = length( find( ins > W_LOWER_BOUND ) );
%     dfWHoldFlag = 0;
%     if mod( it, 10 ) == 0
%         if isempty( find( diff( dfWVec((it-10+1):it))  ~= 0, 1 ) )
%             dfWHoldFlag = 1;
%         end
%     end
    % one of the conditions must be satisfied to stop outer loops
    if abs(LPAry(it+1)-LPAry(it)) <= LP_TOL*LPAry(it) || ...
        ( max( tmp1, tmp2 ) < max( tmp1, tmp2 )*W_TOL &&  tmp3 < D_TOL ) %|| ...
        %( tmp3 < D_TOL && dfWHoldFlag == 1 ) 
        break;
    end
    % the optimization in W or D cannot make the cost function go down.
    if expRec.WStuckFlag == 1 || expRec.DStuckFlag == 1
         break;
    end
end
if ~isempty(filePath)
    fclose(fID);
end
expRec.theta = theta; expRec.lambda = lambda; expRec.phi = phi;
expRec.LPAry = LPAry; expRec.param = param; expRec.aMatrix = aMatrix;
expRec.W = W; expRec.W0 = W0;
expRec.z0 = z0; expRec.z1 = z1; expRec.z2 = z2;
expRec.diffW = tmp1; expRec.diffW0 = tmp2;
expRec.D = D;
expRec.diffD = tmp3;

%expRec.rhoCell = rhoCell; expRec.resRecCell = resRecCell;

end

