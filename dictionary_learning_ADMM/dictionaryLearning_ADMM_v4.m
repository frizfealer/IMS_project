function [ expRec ] = dictionaryLearning_ADMM_v4( inY, initVar, DTemplate, logFY, lambda, theta, phi, aMatrix, BlkDS, filePath, snapFilePath, param, varargin )
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
DHistCell = cell( OUTER_IT_NUM, 1 );


LPAry = zeros( 1, OUTER_IT_NUM+1 );

%rhoCell will consist of MAXIT iteration number of rhoAry
rhoCell = cell( OUTER_IT_NUM, 1 );
%residual Record Cell will consist of MAXIT iteration number of resRecAry
resRecCell = cell(OUTER_IT_NUM, 1);

%% initialize all variable
if isempty(initVar)
    [ outD ] = initD( inY, DTemplate, INIT_METHOD );
    outW = sparse( zeros( mLen, hei*wid ) );
    outW0 = sparse( zeros( hei, wid ) );
else
    outD = initVar.outD;
    outW = initVar.outW;
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

indMap = BlkDS.indMap;

%% managing matlab pool
if matlabpool('size') == 0
    matlabpool( 'open', CLUSTER_NAME, CLUSTER_NUM );
end
%% managing late update for some dictionary elements
% validMap = indMap .* aMatrix;
% yy = inY(:, validMap==1);
[ rSet ] = estimateDRelevant( inY, outD, DTemplate, indMap, LATE_UPDATE_PERCENT, LATE_UPDATE_INTTHRES  );
%% main program start
LPAry(1) = LP_DL_Poiss( aMatrix, inY, outW, outW0, outD, lambda, phi, theta, scaleFactor, logFY );
fprintf( 'parameters: outer iteration number = %d, ', OUTER_IT_NUM );
fprintf( 'ADMM iteration number = %d, Hessian flag for update D = %d, ', ADMM_IT_NUM, HES_FLAG );
fprintf( 'Cluster number used to update = %d, ', CLUSTER_NUM );
fprintf( 'SAVE_TEM_PERIOUD = %d, Inner iteration reduce flag = %d, ', SAVE_TEM_PERIOD, INNER_IT_RED_FLAG );
fprintf( 'LATE_UPDATE_FLAG = %d, LATE_UPDATE_PERCENT=%g, LATE_UPDATE_INTTHRES = %g, ', LATE_UPDATE_FLAG, LATE_UPDATE_PERCENT, LATE_UPDATE_INTTHRES  );
fprintf( ['INIT_D= ', INIT_METHOD, '\n'] );
for it = 1:OUTER_IT_NUM
    fprintf( fID, 'iteration %d\n', it );
    curRhoAry = zeros( ADMM_IT_NUM, 1 );
    curRhoAry(1) = 1;
    resRecAry = zeros( ADMM_IT_NUM, 4 );
    LPAryADMM = zeros( ADMM_IT_NUM, 1 );
    
    % update W, W0, z0, z1
    u0 = sparse( zeros( size( z0 ) ) );
    u1 = sparse( zeros( size( z1 ) ) );
    if INNER_IT_RED_FLAG == 1
        if it < floor(OUTER_IT_NUM*0.8)
            M_ADMM_IT_NUM = round( ADMM_IT_NUM / 2 );
        else
            M_ADMM_IT_NUM = ADMM_IT_NUM;
        end
    else
        M_ADMM_IT_NUM = ADMM_IT_NUM;
    end
    if LATE_UPDATE_FLAG == 1
        if it < floor(OUTER_IT_NUM*0.8)
            relD = outD(:, rSet);
            relW = outW(rSet, :);
            relZ1 = z1(rSet, :);
            relU1 = u1(rSet, :);
            moutD = [relD ones(sLen, 1)];
            tmpEye = eye(length(rSet)+1, length(rSet)+1);
            tmpEye(end, end) = 0;
            CforUW = sparse( [moutD; tmpEye] );
        else
            relD = outD;
            relW = outW;
            relZ1 = z1;
            relU1 = u1;
            tmpEye = eye(mLen+1, mLen+1);
            tmpEye(end, end) = 0;
            moutD = [relD ones(sLen, 1)];
            CforUW = sparse( [moutD; tmpEye] );
        end
    else
        relD = outD;
        relW = outW;
        relZ1 = z1;
        relU1 = u1;
        tmpEye = eye(mLen+1, mLen+1);
        tmpEye(end, end) = 0;
        moutD = [relD ones(sLen, 1)];
        CforUW = sparse( [moutD; tmpEye] );
    end
    prevW = outW; prevW0 = outW0; prevD = outD;
    fprintf( 'M_ADMM_IT_NUM = %d\n', M_ADMM_IT_NUM );
    for itNumADMM = 1:M_ADMM_IT_NUM
        prevWADMM = relW;
        LPAryADMM(itNumADMM) = LP_DL_Poiss( aMatrix, inY, relW, outW0, relD, lambda, phi, theta, scaleFactor, logFY );
        fprintf( 'LP: %g ', LPAryADMM(itNumADMM) );
        tic
        curRho = curRhoAry(itNumADMM);
        fprintf( 'updating W...\t' );
        parfor i = 1:wid*hei
            if  indMap(i)
                [ w_w0 ] = update_w_w0_v2( z0(:,i), u0(:,i), relZ1(:,i), relU1(:,i), CforUW, curRho, [relW(:,i); outW0(i)] );
                relW(:,i) = w_w0(1:end-1); outW0(i) = w_w0(end);
            end
        end
        fprintf( 'updating z0...\t' );
        prevZ0 = z0;    %for computing the dual residual
        parfor i = 1:hei*wid
            if indMap(i) && aMatrix(i)
                [ z0(:, i) ] = updatez0_v3( aMatrix(i), inY(:,i), z0(:,i), relD, relW(:, i), outW0(i), u0(:,i), curRho, scaleFactor );
            end
        end
        
        fprintf( 'updating z1...\n' );
%         prevZ1 = z1;    %for computing the dual residual
        for j = 1:BlkDS.blkNum
            loc = BlkDS.B2GMap{j};
            curW = relW(:, loc);
            curU1 = relU1(:, loc);
            tmpZ1 = zeros( size(relW, 1), length( loc ) );
            curR = Rblk{j};
            parfor i = 1:size(relW, 1)
%                 [ tmpZm ] = updatez1_m( [], [], curRho, curW(i, :), curU1(i, :), lambda+1e-8, theta+1e-8, curR, []);
                tmpZm = updatez1_m_ridgePen_SLEP( curRho, curW(i, :), curU1(i, :), lambda, theta,curR );
                tmpZ1(i, :) = tmpZm;
            end
            relZ1(:, loc) = tmpZ1;
        end
    %     after = LP_DL_Poiss( aMatrix, inY, outW, outW0, outD, lambda, 0, theta );
    %     if abs(after-before) > 1e-6 && after > before
    %         keyboard();
    %     end
    %     before = after;

        preY = relD*relW(:,:) + repmat( outW0(:)', sLen, 1);
        %primal residual r1
        r1 = z0 - preY;
%         [val, o] = max(r1(:));
%         fprintf( '%g %g %g\n', val,  o, full(z0(o))); 
%         for i = 1:wid*hei
%             if ~indMap(i) && ~isempty( find( r1(:,i) ~= 0, 1 ) )
%                 fprintf( 'loc: %d', i );
%                 r1(:, i) = 0;
%             end
%         end
        u0 = u0 + curRho*r1;
        %compute the primal residual all
        rf = r1;
        r2 = relZ1- relW;
%         for i = 1:wid*hei
%             if ~indMap(i) && ~isempty( find( r2(:,i) ~= 0, 1 ) )
%                fprintf( 'loc: %d', i );
%                 r2(:, i) = 0;
%             end
%         end
        
        %                 fprintf( 'r2 < 0 : %d\n', length(find(r2b(:)<0)) );
        relU1 = relU1 + curRho* r2;
        rf = [rf; r2];
        rfn2 = norm(rf(:), 2);
        
        %compute the dual variable to see if it is converge
%         dualVar = [u0; u1{1}];
%         dualVar = norm(dualVar(:), 2);
        
        %dual residual s0
        s0 = curRho*( relD'*( z0 - prevZ0 ) );
%         for i = 1:wid*hei
%             if ~indMap(i) && ~isempty( find( s0(:,i) ~= 0, 1 ) )
%                 fprintf( 'loc: %d', i );
%                 s0(:,i) = 0;
%             end
%         end
        s0 = [s0; sum( z0 -prevZ0 )];
        s0n2 = norm( s0(:), 2 );
        
        %compute the more dual residual values needed at the block condtion
%         duResMore = [];
%         for i = 1:hei*wid
%             bID = BlkDS.mapW2B{i};
%             if length(bID) > 1
%                 tmp = 0;
%                 for j = 2:length(bID)
%                     tmp = tmp + curRho*( ( z1{bID(j)}(:, i)- prevZ1{bID(j)}(:, i) ) );
%                 end
%                 duResMore = [duResMore; tmp];
%             end
%         end
%         duResMore = norm(duResMore, 2);

        %converge condition
        %primal residual condition
        EPS_ABS = 1e-4;
        EPS_REL = 1e-3;
        
        epsPri = sqrt( size(rf, 1)*size(rf, 2) )*EPS_ABS;
        relEpsPri = max( [norm( z0(:), 2 ) norm( preY(:) ) norm( relZ1(:), 2 ), norm( relW(:), 2 )] );
        epsPri = epsPri + relEpsPri*EPS_REL;
        
        %dual residual condition
%         tmp = outD'*u0;
%         for i = 1:BlkDS.blkNum
%             if dataMask(i) == 0
%                 continue;
%             else
%                 tmp = tmp + u1{i};
%             end
%         end
%         tmp = [tmp; sum(u0) ];
%         epsDual = sqrt( size(s0, 1) * size(s0, 2) )*1e-4 + 1e-3 * norm( tmp(:), 2);
        tmp = relD'*u0;
        tmp = [tmp; sum(u0)];
        maxRel = max( norm( tmp(:), 2 ), norm( relU1(:), 2 ) );
        epsDual = sqrt( size(s0, 1) * size(s0, 2) )*EPS_ABS + maxRel*EPS_REL;
        
        resRecAry(itNumADMM, :) = [rfn2, epsPri, s0n2, epsDual];
        tmp = max( abs( relW(:) -  prevWADMM(:) ) );
        fprintf( '%d: %g %g %g %g %g %g\n',itNumADMM, rfn2, epsPri, s0n2, epsDual, curRhoAry(itNumADMM), full(tmp) );
        if ( rfn2 < epsPri && s0n2 < epsDual ) ||...
                ( max( abs(rf(:)) ) < 1e-3 && max( abs(s0(:)) ) < 1e-3 ) || ...
                ( tmp < 1e-3 )
            break;
        end

        %update the next stage rho
%         if rfn2 > 2*s0n2
%             curRhoAry(itNumADMM+1) = 1.5 * curRhoAry(itNumADMM);
        if rfn2 > 10*s0n2
            curRhoAry(itNumADMM+1) = 2 * curRhoAry(itNumADMM);
%         elseif s0n2 > 2*rfn2
%             curRhoAry(itNumADMM+1) = 0.677 * curRhoAry(itNumADMM);
        elseif s0n2 > 10*rfn2
            curRhoAry(itNumADMM+1) = 0.5 * curRhoAry(itNumADMM);
        else
            curRhoAry(itNumADMM+1) = curRhoAry(itNumADMM); 
        end
        toc
    end
    rhoCell{it} = curRhoAry;
    resRecCell{it} = resRecAry;

    
%     [B, FitInfo] = lassoglm( [outD], inY, 'poisson', 'Lambda', 0, 'RelTol', 1e-8 ); 
%     cvx_begin
%         variable uW(mLen+1);
%         maximize( sum( inY .*([outD ones(sLen, 1)]*uW) -  exp([outD ones(sLen, 1)]*uW) ) - lambda*norm(uW(1:end-1), 1) );
%         subject to
%             uW >= 0;
%     cvx_end
 
    %update D
    if INNER_IT_RED_FLAG == 1
        if it < floor(ADMM_IT_NUM*0.8)
            M_UP_D_IT_NUM = round( UP_D_IT_NUM / 2 );
        else
            M_UP_D_IT_NUM = UP_D_IT_NUM;
        end
    else
        M_UP_D_IT_NUM = UP_D_IT_NUM;
    end
    validMap = indMap .* aMatrix;
    if isempty( deInfo )
        if LATE_UPDATE_FLAG == 1 && it < floor(OUTER_IT_NUM*0.8)
           %[ relD ] = updateD_v5( inY, relW, outW0, relD, DTemplate(:, rSet), validMap, HES_FLAG, phi, scaleFactor, M_UP_D_IT_NUM );
            [ relD ] = updateD_v8_ipopt( inY, relW, outW0, relD, DTemplate(:, rSet), validMap, HES_FLAG, phi, scaleFactor, M_UP_D_IT_NUM );
        else
%            [ relD ] = updateD_v5( inY, relW, outW0, relD, DTemplate, validMap, HES_FLAG, phi, scaleFactor, M_UP_D_IT_NUM );
            [ relD ] = updateD_v8_ipopt( inY, relW, outW0, relD, DTemplate, validMap, HES_FLAG, phi, scaleFactor, M_UP_D_IT_NUM );
        end
    end
%     
%     z0 = sparse( z0Init);
%     z1 = cell(BlkDS.blkNum, 1);
%     for i = 1:length(z1)
%         z1{i} = zeros( mLen, hei*wid );
%     end
%     u0 = sparse(zeros(sLen, hei*wid));
%     u1 = cell(BlkDS.blkNum, 1);
%     for i = 1:length(u1)
%         u1{i} = sparse( mLen, hei*wid );
%     end
    if LATE_UPDATE_FLAG == 1
        if it < floor(OUTER_IT_NUM*0.8)
            outD(:, rSet) = relD;
            outW(rSet, :) = relW;
            z1(rSet, :) = relZ1;
            u1(rSet, :) = relU1;
        else
            outD = relD;
            outW = relW;
            z1 = relZ1;
            u1 = relU1;
        end
    else
        outD = relD;
        outW = relW;
        z1 = relZ1;
        u1 = relU1;
    end
    
    LPAry(it+1) = LP_DL_Poiss( aMatrix, inY, outW, outW0, outD, lambda, phi, theta, scaleFactor, logFY );
    tmp1 = max( abs( outW(:)-prevW(:) ) );
    tmp2 = max( abs( outW0(:)-prevW0(:) ) );
    tmp3 = max( abs( outD(:)-prevD(:) ) );
    if mod(it, SAVE_TEM_PERIOD) == 0
        expRec.outD = outD; expRec.outW = outW; expRec.outW0 = outW0;
        expRec.LPAry = LPAry; expRec.z0 = z0; expRec.z1 = z1; 
        expRec.u0 = u0; expRec.u1 = u1;
        expRec.LPAryADMM = LPAryADMM;
        expRec.theta = theta; expRec.lambda = lambda; expRec.phi = phi;
        expRec.diffW = tmp1; expRec.diffW0 = tmp2; expRec.diffD = tmp3;
%         expRec.rhoCell = rhoCell; expRec.resRecCell = resRecCell;
        save( snapFilePath, 'expRec' );
    end
    DHistCell{it} = outD;
    if ~isempty( D_HIST_PATH )
        save( D_HIST_PATH, 'DHistCell' );
    end
    if abs(LPAry(it+1)-LPAry(it)) <= 1e-6 || ...
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

