function [ expRec ] = dictionaryLearning_ADMM_v5( inY, initVar, DTemplate, logFY, lambda, theta, phi, aMatrix, BlkDS, filePath, snapFilePath, param, varargin )
%dictionaryLearning_ADMM_4 dictionary learning with ADMM formulation
%with set z0 = D*W+W0, z1 = w1, z2 = Ez1;
%Dictioanry Learning with Poisson distribution
%input:
% inY [s h w]
% initVar data structure with all initialization, has the fields:
% initVar.outW initVar.outW0 initVar.outD 
% initVar.u0 initVar.u1 initVar.z0 initVar.z1
%DTemplate the dictionary pattern you want to have [s m]
%We choose the appropriate dictionary pattern before this function
%lambda, theta: scalar
%filePath: the output file path, if output to screen, set []
%aMatrix [h w] a Matrix indicate which grid is held out, 1 means trainning,
%0 means testing
%snapFilePath: file path for snapshot every 20 iteraionts
%dataMask: the region we donot want not need to update w, z1, and u1, z0,
%u0: with size [blkNum 1], where 1 means need counts and 0 the opposite
%BlkDS, structures with two fields,
%   bIdx is the block to location index
%   mapW2B is the location to block index
%   blkNum the number of the blocks
%last variable can be iteration number
%output:
%expRec is a structure contains 
%LPAry [itNum+1, 1], a log posterior record after each iteration, the first
%one is the log posterior from orignal iteration
%outW, outW0, outD, z0, z1, u0, u1
%rho, theta, lambda
[sLen, hei, wid] = size(inY);
[mLen] = size(DTemplate ,2);

%% initializing z0 as log(inY)
z0 = log(inY);
z0(z0==-inf)=0;
%if not in traing set, we shoud not initialize z0 according to it.
%take the average values
tmp = sum( z0(:, :), 2 ) / ( length( find( BlkDS.indMap == 1 ) ) );
for i = 1:hei*wid
    %if is in data area, but not in training set, 
    %set it to average spectrum
    if BlkDS.indMap(i) && ~aMatrix(i)
        z0(:, i) = tmp;
    end
end
scaleFactor =  1 / ( max( inY(:) ) / max( z0(:) ) ) * 100;
% scaleFactor = 1;
%% change to sparse matrix type, 2D, initialize z0, z1, z2, u0, u1, u2
z0 = sparse(z0(:, :));
z1 = sparse( mLen, hei*wid );
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
    ADMM_IT_NUM = 100;
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
    CLUSTER_NUM = 12;
end
if isfield( param, 'SAVE_TEM_PERIOD' )
    SAVE_TEM_PERIOD = param.SAVE_TEM_PERIOD;
else
    SAVE_TEM_PERIOD = 5;
end
if isfield( param, 'INNER_IT_RED_FLAG' )
    INNER_IT_RED_FLAG = param.INNER_IT_RED_FLAG;
else
    INNER_IT_RED_FLAG = 0;
end
if isfield( param, 'LATE_UPDATE_FLAG' )
    LATE_UPDATE_FLAG = param.LATE_UPDATE_FLAG;
else
    LATE_UPDATE_FLAG = 1;
end
if isfield( param, 'LATE_UPDATE_PERCENT' )
    LATE_UPDATE_PERCENT = param.LATE_UPDATE_PERCENT;
else
    LATE_UPDATE_PERCENT = 0.2;
end
if isfield( param, 'LATE_UPDATE_INTTHRES' )
    LATE_UPDATE_INTTHRES = param.LATE_UPDATE_INTTHRES;
else
    LATE_UPDATE_INTTHRES = 0.8;
end
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

LPAry = zeros( 1, OUTER_IT_NUM+1 );

%rhoCell will consist of MAXIT iteration number of rhoAry
rhoCell = cell( OUTER_IT_NUM, 1 );
%residual Record Cell will consist of MAXIT iteration number of resRecAry
resRecCell = cell(OUTER_IT_NUM, 1);

%% initialize all variable
if isempty(initVar)
    [ outD ] = initD( inY, DTemplate, INIT_METHOD, D_ION_NAME );
    DHistCell{1} = outD;
    outW = sparse( zeros( mLen, hei*wid ) );
    WHistCell{1} = outW;
    outW0 = sparse( zeros( hei, wid ) );
else
    outD = initVar.outD;
    DHistCell{1} = outD;
    outW = initVar.outW;
    WHistCell{1} = outW;
    outW0 = initVar.outW0;
    z0 = initVar.z0; z1 = initVar.z1;
end

%% for debuging
if length(varargin)>=1
    deInfo = varargin{1};
    outD = deInfo.outD;
else
    deInfo = [];
end

%% contructing genSparseGroupingMatrix
[Rall] = genSparseGroupingMatrix( hei, wid, 1 );
for i = 1:BlkDS.blkNum
    Rblk{i} = [];
    Rb = Rall(:, BlkDS.B2GMap{i});
    ins = sum(abs(Rb), 2);
    Rb(ins<2,:) = [];
    Rblk{i} = Rb;
end

%% managing matlab pool
if matlabpool('size') == 0
    matlabpool( 'open', CLUSTER_NAME, CLUSTER_NUM );
end
%% managing late update for some dictionary elements
% validMap = indMap .* aMatrix;
% yy = inY(:, validMap==1);
if LATE_UPDATE_FLAG == 1
    [ rSet ] = estimateDRelevant( inY, outD, DTemplate, BlkDS.indMap, LATE_UPDATE_PERCENT, LATE_UPDATE_INTTHRES, D_ION_NAME  );
end
%% main program start
MEAN_FLAG = 0;
LPAry(1) = LP_DL_Poiss( aMatrix, inY, outW, outW0, outD, lambda, phi, theta, scaleFactor, logFY, MEAN_FLAG );
fprintf( 'parameters: outer iteration number = %d, ', OUTER_IT_NUM );
fprintf( 'ADMM iteration number = %d, D-update iteration number = %d, Hessian flag for update D = %d, ', ADMM_IT_NUM, UP_D_IT_NUM, HES_FLAG );
fprintf( 'Cluster number used to update = %d, ', CLUSTER_NUM );
fprintf( 'SAVE_TEM_PERIOUD = %d, Inner iteration reduce flag = %d, ', SAVE_TEM_PERIOD, INNER_IT_RED_FLAG );
fprintf( 'LATE_UPDATE_FLAG = %d, LATE_UPDATE_PERCENT=%g, LATE_UPDATE_INTTHRES = %g, ', LATE_UPDATE_FLAG, LATE_UPDATE_PERCENT, LATE_UPDATE_INTTHRES  );
fprintf( ['INIT_D= ', INIT_METHOD, '\n'] );
for it = 1:OUTER_IT_NUM
    fprintf( fID, 'iteration %d\n', it );
    if INNER_IT_RED_FLAG == 1
        if it < floor(OUTER_IT_NUM*0.8)
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
    if LATE_UPDATE_FLAG == 1
        if it < floor(OUTER_IT_NUM*0.8)
            relD = outD(:, rSet);
            relW = outW(rSet, :);
            relZ1 = z1(rSet, :);
        else
            relD = outD;
            relW = outW;
            relZ1 = z1;
        end
    else
        relD = outD;
        relW = outW;
        relZ1 = z1;
    end
    prevW = outW; prevW0 = outW0; prevD = outD;
    %update W
    [WResStruct] = updateW_ADMM( inY, relD, relW, outW0, z0, relZ1, aMatrix, BlkDS, M_ADMM_IT_NUM, lambda, phi, theta, scaleFactor, logFY, Rblk );
    
    relW = WResStruct.W;
    outW0 = WResStruct.W0;
    z0 = WResStruct.z0;
    relZ1 = WResStruct.z1;
    LPAryADMM = WResStruct.LPAry;
    rhoCell{it} = WResStruct.rhoAry;
    resRecCell{it} = WResStruct.resAry;
    
    %% if W is too low, make corresponding D to zero
    if LATE_UPDATE_FLAG == 1 && it < floor(OUTER_IT_NUM*0.8)
        relDTemplate = DTemplate(:, rSet);
    else
        relDTemplate = DTemplate;
    end
    W_lOWERBOUND = 1e-4;
    for i = 1:size(relW, 1)
        if max(relW(i,:)) < W_lOWERBOUND
            relD(relDTemplate(:,i)==1, i) = 1e-6;
        end
    end

    validMap = BlkDS.indMap .* aMatrix;
    if isempty( deInfo )
%         if LATE_UPDATE_FLAG == 1 && it < floor(OUTER_IT_NUM*0.8)
%            %[ relD ] = updateD_v5( inY, relW, outW0, relD, DTemplate(:, rSet), validMap, HES_FLAG, phi, scaleFactor, M_UP_D_IT_NUM );
%             [ relD ] = updateD_v8_ipopt( inY, relW, outW0, relD, DTemplate(:, rSet), validMap, HES_FLAG, phi, scaleFactor, M_UP_D_IT_NUM );
%         else
% %            [ relD ] = updateD_v5( inY, relW, outW0, relD, DTemplate, validMap, HES_FLAG, phi, scaleFactor, M_UP_D_IT_NUM );
%             [ relD ] = updateD_v8_ipopt( inY, relW, outW0, relD, DTemplate, validMap, HES_FLAG, phi, scaleFactor, M_UP_D_IT_NUM );
%         end
        [ relD ] = updateD_v8_ipopt( inY, relW, outW0, relD, relDTemplate, validMap, HES_FLAG, phi, scaleFactor, M_UP_D_IT_NUM );
    end

    if LATE_UPDATE_FLAG == 1
        if it < floor(OUTER_IT_NUM*0.8)
            outD(:, rSet) = relD;
            outW(rSet, :) = relW;
            z1(rSet, :) = relZ1;
        else
            outD = relD;
            outW = relW;
            z1 = relZ1;
        end
    else
        outD = relD;
        outW = relW;
        z1 = relZ1;
    end
    
    LPAry(it+1) = LP_DL_Poiss( aMatrix, inY, outW, outW0, outD, lambda, phi, theta, scaleFactor, logFY, MEAN_FLAG );
    tmp1 = max( abs( outW(:)-prevW(:) ) );
    tmp2 = max( abs( outW0(:)-prevW0(:) ) );
    tmp3 = max( abs( outD(:)-prevD(:) ) );
    if mod(it, SAVE_TEM_PERIOD) == 0
        expRec.outD = outD; expRec.outW = outW; expRec.outW0 = outW0;
        expRec.LPAry = LPAry; expRec.z0 = z0; expRec.z1 = z1; 
        expRec.LPAryADMM = LPAryADMM;
        expRec.theta = theta; expRec.lambda = lambda; expRec.phi = phi;
        expRec.diffW = tmp1; expRec.diffW0 = tmp2; expRec.diffD = tmp3;
        expRec.param = param;
%         expRec.rhoCell = rhoCell; expRec.resRecCell = resRecCell;
        save( snapFilePath, 'expRec' );
    end
    if ~isempty( D_HIST_PATH )
        DHistCell{it+1} = outD;
        save( D_HIST_PATH, 'DHistCell' );
    end
    if ~isempty( W_HIST_PATH )
        WHistCell{it+1} = outW;
        save( W_HIST_PATH, 'WHistCell', '-v7.3' );
    end
    % ALL conditions must be satisfied to stop outer loops
    if abs(LPAry(it+1)-LPAry(it)) <= 1e-6 && ...
        ( max( tmp1, tmp2 ) < 1e-3 &&  tmp3 < 1e-3 )
        break;
    end
end
if ~isempty(filePath)
    fclose(fID);
end
expRec.outD = outD; expRec.outW = outW; expRec.outW0 = outW0;
expRec.LPAry = LPAry; expRec.z0 = z0; expRec.z1 = z1;
expRec.LPAryADMM = LPAryADMM;
expRec.theta = theta; expRec.lambda = lambda; expRec.phi = phi; 
expRec.diffW = tmp1; expRec.diffW0 = tmp2; expRec.diffD = tmp3;
%expRec.rhoCell = rhoCell; expRec.resRecCell = resRecCell;

end

