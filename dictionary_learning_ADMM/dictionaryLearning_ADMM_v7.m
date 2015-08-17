function [ DL ] = dictionaryLearning_ADMM_v7( InputData, DL, verboseW, verboseD, snapFilePath )
%--------------------------------------------------------------------------
% dictionaryLearning_ADMM_v7: dictionary learning with ADMM formulation 
% version 7 with dictionary learning data structure
%--------------------------------------------------------------------------
% DESCRIPTION:
%   The main function of the algorithm. Poisson distribution.
%   ADMM setting: z0 = D*W+W0, z1 = w1, z2 = RW';
% INPUT ARGUMENTS:
%   InputData, the InputData data structure.
%   DL, the dictionary learning data structure.
%   verboseW, a flag indicate whether to show the messages of updating W.
%   verboseD, a flag indicate whether to show the messages of updating D.
%   snapFilePath, the output file path for the snapshot of the model.
% OUTPUT ARGUMENTS:
%   DL, the updated dictionary leanring data structure.

%% pour out data
Y = InputData.dataCube;
aMatrix = InputData.aMatrix;
scaleFactor = InputData.scaleFactor;
logFY = InputData.logFY;

%% set all parameters
% if isfield( param, 'CLUSTER_NUM' ) %deprecated
%     CLUSTER_NUM = param.CLUSTER_NUM;
% else
%     CLUSTER_NUM = 1;
% end

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
% if isfield( param, 'LATE_UPDATE_IT_PREC' ) %deprecated
%     LATE_UPDATE_IT_PREC = param.LATE_UPDATE_IT_PREC;
% else
%     LATE_UPDATE_IT_PREC = 0.8;
% end
%for second-round dictionary learning, under construction
% if isfield( param, 'newWInfo' )
%     newWInfo = param.newWInfo;
% else
%     newWInfo = [];
% end
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


DHistCell = cell( DL.itNum+1, 1 );
WHistCell = cell( DL.itNum+1, 1 );

%rhoCell will consist of MAXIT iteration number of rhoAry
% rhoCell = cell( OUTER_IT_NUM, 1 );
%residual Record Cell will consist of MAXIT iteration number of resRecAry
% resRecCell = cell(OUTER_IT_NUM, 1);
% kappa = 1e-4;
%% managing matlab pool
% if matlabpool('size') == 0 && CLUSTER_NUM > 1
%     matlabpool( 'open', CLUSTER_NAME, CLUSTER_NUM );
% end
%%managing late update for some dictionary elements
% if LATE_UPDATE_FLAG == 1
%     [ rSet ] = estimateDRelevant( inY, outD, DTemplate, BlkDS.indMap, LATE_UPDATE_PERCENT, LATE_UPDATE_INTTHRES, D_ION_NAME  );
% end
%% main program start
DL.LPAry(1) = LP_DL_Poiss( DL.W.LINK_FUNC, aMatrix, Y, DL.W.W, DL.W.W0, DL.D.D, DL.W.lambda, DL.D.phi, DL.W.theta, scaleFactor, logFY, 0 );
fprintf( 'parameters: outer iteration number = %d \n', DL.itNum );
fprintf( 'ADMM iteration number = %d, ...D-update iteration number = %d. SAVE_TEM_PERIOUD = %d.\n', DL.W.itNum, DL.D.itNum, DL.SAVE_TEM_PERIOD );
% fprintf( 'Cluster number used to update = %d, ', CLUSTER_NUM );
fprintf( ['INIT_D= ', DL.D.INIT_METHOD, '.\n'] );
fprintf( ['LINK_FUNC= ', DL.W.LINK_FUNC, '.\n'] );

fprintf( 'LP: %g\n', DL.LPAry(1) );
for it = 1:DL.itNum
    fprintf( 'iteration %d\n', it );
%     if INNER_IT_RED_FLAG == 1
%         if it < floor(OUTER_IT_NUM*LATE_UPDATE_IT_PREC)
%             M_ADMM_IT_NUM = round( ADMM_IT_NUM / 2 );
%             M_UP_D_IT_NUM = round( UP_D_IT_NUM / 2 );
%         else
%             M_ADMM_IT_NUM = ADMM_IT_NUM;
%             M_UP_D_IT_NUM = UP_D_IT_NUM;
%         end
%     else
%         M_ADMM_IT_NUM = ADMM_IT_NUM;
%         M_UP_D_IT_NUM = UP_D_IT_NUM;
%     end
    prevW = DL.W.W; prevW0 = DL.W.W0; prevD = DL.D.D;
    %% update W
    [DL.W] = updateW_ADMM_v7( InputData, DL.D.D, DL.W, verboseW );
    if DL.W.stuckFlag == 1
        break;
    end
    
    %% update D
    before = LP_DL_Poiss( DL.W.LINK_FUNC, aMatrix, Y, DL.W.W, DL.W.W0, DL.D.D, DL.W.lambda, DL.D.phi, DL.W.theta, scaleFactor, logFY, 0 );
    [ DL.D, ~ ] = updateD_v9_ipopt( InputData, DL.W.W, DL.W.W0, DL.D, verboseD );
    after = LP_DL_Poiss( DL.W.LINK_FUNC, aMatrix, Y, DL.W.W, DL.W.W0, DL.D.D, DL.W.lambda, DL.D.phi, DL.W.theta, scaleFactor, logFY, 0 );
    
    cnt = 0;
    while after-before > 0
        if cnt == 5
            DL.D.stuckFlag = 1;
            break;
        end
        cnt = cnt + 1;
        [ DL.D, ~ ] = updateD_v9_ipopt( InputData, DL.W.W, DL.W.W0, DL.D, verboseD );
        after = LP_DL_Poiss( DL.W.LINK_FUNC, aMatrix, Y, DL.W.W, DL.W.W0, DL.D.D, DL.W.lambda, DL.D.phi, DL.W.theta, scaleFactor, logFY, 0 );
    end
    if DL.D.stuckFlag == 1
        break;
    end
    
    %% compute the termination criteria
    DL.LPAry(it+1) = LP_DL_Poiss( DL.W.LINK_FUNC, aMatrix, Y, DL.W.W, DL.W.W0, DL.D.D, DL.W.lambda, DL.D.phi, DL.W.theta, scaleFactor, logFY, 0 );
    tmp1 = max( abs( DL.W.W(:)-prevW(:) ) );
    tmp2 = max( abs( DL.W.W0(:)-prevW0(:) ) );
    tmp3 = max( abs( DL.D.D(:)-prevD(:) ) );
    if DL.SAVE_TEM_PERIOD ~= -1 && mod(it, DL.SAVE_TEM_PERIOD) == 0
        save( snapFilePath, 'DL' );
    end
    if DL.D_HIST_FLAG == 1
        DHistCell{it+1} = DL.D.D;
    end
    if DL.W_HIST_FLAG == 1
        WHistCell{it+1} = DL.W.W;
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
    if abs(DL.LPAry(it+1)-DL.LPAry(it)) <= DL.LP_TOL*DL.LPAry(it) || ...
        ( max( tmp1, tmp2 ) < max( [DL.W.W(:); DL.W.W0(:)] )*DL.W_TOL ||...
        tmp3 < DL.D_TOL ) %|| ...
        %( tmp3 < D_TOL && dfWHoldFlag == 1 ) 
        break;
    end
    fprintf( 'LP: %g diff LP: %g, tol LP: %g, diff W: %g, tol W: %g, diff D: %g, tol D: %g\n', ...
        DL.LPAry(it+1), ...
        abs(DL.LPAry(it+1)-DL.LPAry(it)), DL.LP_TOL*DL.LPAry(it), ...
       max( tmp1, tmp2 ), max( [DL.W.W(:); DL.W.W0(:)] )*DL.W_TOL, ...
       tmp3, DL.D_TOL );
end
DL.DHistCell = DHistCell;
DL.WHistCell = WHistCell;
end

