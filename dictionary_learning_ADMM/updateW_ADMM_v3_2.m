function [WResStruct] = updateW_ADMM_v3_2( LINK_FUNC, Y, D, aMatrix, itNum, lambda, theta, L1Flag, logFY, initVar, scaleFactor, varargin )

%% construct metaData
BlkDS = conBLKDS( Y );
[sLen, hei, wid] = size( Y );
% nLen = hei*wid;
[Rall] = genSparseGroupingMatrix( hei, wid, 1 );
for i = 1:BlkDS.blkNum
    Rblk{i} = [];
    Rb = Rall(:, BlkDS.B2GMap{i});
    ins = sum(abs(Rb), 2);
    Rb(ins<2,:) = [];
    Rblk{i} = Rb;
end

%% variable setting
newWInfo = [];
if ~isempty(varargin)  
    if length(varargin) >= 3
        WHistFlag = varargin{1};
        Wtol = varargin{2};
        D_LOWER_BOUND = varargin{3};
    else
        Wtol = 5*1e-3;
        D_LOWER_BOUND = 1e-2;
    end
    if length(varargin) >= 4
        newWInfo = varargin{4};
    end
    if length(varargin) >= 5
        kappa = varargin{5};
    end
else
    Wtol = 5*1e-3;
    D_LOWER_BOUND = 1e-2;
end

fprintf( 'iteration number for W_update_ADMM = %d\n', itNum );
temp = (aMatrix==0)&(BlkDS.indMap==1);
fprintf( 'held-out sample number = %d, training sample number = %d\n', length( find( temp == 1 ) ), length( find( aMatrix == 1 ) ) );
fprintf( 'lambda = %d, theta = %d\n', lambda, theta );
%% modifying D
mLen = size(D, 2);
ins = max( D, [], 1 );
rMIdx = find( ins > D_LOWER_BOUND );
D = D(:, rMIdx);
oMLen = mLen;
mLen = length( rMIdx );
if ~isempty( newWInfo )
    for i = 1:length(newWInfo)
        if BlkDS.indMap(i) == 1
            newWInfo{i} = intersect( newWInfo{i}, rMIdx );
        end
    end
end
%% initialize all variables( W, W0, z0, z1, z2, u0, u1, u2)
if isempty( initVar )
    %setting W, W0, z0, z1, z2
    if strcmp( LINK_FUNC, 'log' ) == 1 
        z0 = log(Y); z0(z0==-inf)=0; %initialize zo as log(Y)
    elseif strcmp( LINK_FUNC, 'identity' ) == 1
        z0 = Y;
        z0( z0 == 0 ) = 1e-8;
    elseif strcmp( LINK_FUNC, 'log_gaussain' ) == 1
        z0 = log(Y); z0(z0==-inf)=0; %initialize zo as log(Y)    
    elseif strcmp( LINK_FUNC, 'negative_binomial' ) == 1
        z0 = log(Y); z0(z0==-inf)=0;
    end
    %if not in traing set, we shoud not initialize z0 according to it.
    %take the average values
    idx = intersect( find(aMatrix==1), find(BlkDS.indMap==1) );
    tmp = sum( z0(:, idx), 2 ) / ( length( idx ) );
    for i = 1:hei*wid
        %if is in data area, but not in training set,
        %set it to average spectrum
        if BlkDS.indMap(i) == 1 && aMatrix(i) == 0
            z0(:, i) = tmp;
        end
    end
    z1 = zeros( mLen, wid*hei ); 
    z2 = cell( BlkDS.blkNum, 1 );
    for i = 1:BlkDS.blkNum
        z2{i} = zeros( length(Rblk{i}), mLen );
    end
    W = zeros( mLen, wid*hei ); W0 = zeros( hei, wid );
    %construct a initVar variable, the size is the same as the input size
    initVar.z0 = z0; initVar.z1 = zeros( oMLen, wid*hei );
    for i = 1:BlkDS.blkNum
        initVar.z2{i} = zeros( length(Rblk{i}), oMLen );
    end
    initVar.W = zeros( oMLen, wid*hei ); initVar.W0 = zeros( hei, wid );
else
    z0 = initVar.z0(:, :);
    if strcmp( LINK_FUNC, 'identity' ) == 1
        z0( z0 == 0 ) = 1e-8;
    end
    z1 = initVar.z1(rMIdx, :); 
    W = initVar.W(rMIdx, :); W0 = initVar.W0;
    z2 = cell( BlkDS.blkNum, 1 );
    for i = 1:BlkDS.blkNum
        z2{i} = initVar.z2{i}(:, rMIdx);
    end
end
%setting u0, u1, u2 and init.u0, u1. u2
u0 = sparse( zeros( size( z0(:, :) ) ) ); initVar.u0 = u0;
u1 = sparse( zeros( size( z1(:, :) ) ) ); initVar.u1 = sparse( zeros( oMLen, wid*hei ) );
u2 = cell( BlkDS.blkNum, 1 ); initVar.u2 = cell( BlkDS.blkNum, 1 );
for i = 1:BlkDS.blkNum
    u2{i} = sparse( zeros( size( z2{i}(:, :) ) ) );
    initVar.u2{i} = sparse( zeros( length(Rblk{i}), oMLen ) );
end

%% initialize others
LPAryADMM = zeros( itNum+1, BlkDS.blkNum );
resRecAry = zeros( itNum, 4, BlkDS.blkNum );
curRhoAry = zeros( itNum, BlkDS.blkNum );
curRhoAry(1, :) = 1;

moutD = [D ones(sLen, 1)];
tmpEye = eye(mLen+1, mLen+1);
tmpEye(end, end) = 0;
aPartofCforUW = sparse( [moutD; tmpEye] );
stuckFlag = 0;
oITNum = itNum;
%% main optimization process
for j = 1:BlkDS.blkNum %for a block
    LPAryADMM(1, j) = LP_DL_Poiss( LINK_FUNC, aMatrix, Y, W, W0, D, lambda, 0, theta, scaleFactor, logFY, 0 );
    %initialize variables in a block
    loc = BlkDS.B2GMap{j};
    curY = Y(:, loc);
    curW = W(:, loc);
    curW0 = W0(loc);
    curZ0 = z0(:, loc);
    curU0 = u0(:, loc);
    curZ1 = z1(:, loc);
    curU1 = u1(:, loc);
    curZ2 = z2{j};
    curU2 = u2{j};
    cRblk = Rblk{j};
    CforUW = replicateC( aPartofCforUW, length(loc) );
    nrPart1 = size(CforUW, 1);
    RforUW = genSparseGroupingMatrix2( cRblk, mLen, 1 );
    CforUW = rowCombSparseMatrix( CforUW, RforUW );
    if ~isempty(newWInfo)
        [CforUW, target] =  modifyCWithSparsity( CforUW, newWInfo, loc, mLen+1, nrPart1 );
    end
    fprintf( 'LP %g\n', LPAryADMM(1, j) );
    itNumADMM = 1;
    while itNumADMM <= itNum
%         tic;
        preW = curW;
        %check to see at least LPAryADMM need to larger than the origanl
        %LPAryADMM LPAryADMM(1, j) 
        curRho = curRhoAry(itNumADMM, j);
        %% update W
        fprintf( 'updating W... ' );
        [ w_w0 ] = update_w_w0_v4( curZ0, curU0, curZ1, curU1, curZ2, curU2, CforUW, curRho, [curW; curW0(:)'], 1 );
        tmp = reshape(w_w0, mLen+1, length(loc));
        curW = tmp(1:(end-1), :);
        curW0 = tmp(end,:);
        %% update Z0
        fprintf( 'updating z0... ' );
        locAlpha0 = intersect( loc, find( aMatrix(:) == 0 ) );
        locAlpha0 = find( ismember( loc, locAlpha0 ) );
        locAlpha1 = intersect( loc, find( aMatrix(:) == 1 ) );
        locAlpha1 = find( ismember( loc, locAlpha1 ) );
        %update for held-out data
        preZ0 = curZ0;
        [ curZ0(:, locAlpha0) ] = updatez0_v4( 0, LINK_FUNC, curY(:, locAlpha0), curZ0(:, locAlpha0), D, ...
            curW(:, locAlpha0), curW0(locAlpha0), curU0(:, locAlpha0), curRho, scaleFactor );
        %update for training data
        [ curZ0(:, locAlpha1) ] = updatez0_v4( 1, LINK_FUNC, curY(:, locAlpha1), curZ0(:, locAlpha1), D, ...
            curW(:, locAlpha1), curW0(locAlpha1), curU0(:, locAlpha1), curRho, scaleFactor, [], [] );
        %% update z1
        fprintf( 'updating z1... ' );
        preZ1 = curZ1;
        [ curZ1 ] = updatez1_v2( curZ1, curU1, curW, curRho, lambda );
        %% update z2
        fprintf( 'updating z2... ' );
        preZ2 = curZ2;
        [ curZ2 ] = updatez2_v1( curZ2, curU2, curRho, cRblk, curW, theta, L1Flag );
        %% update u0, u1
        fprintf( 'updating u0, u1, u2...\n ' );
        preY = D*curW(:,:) + repmat( curW0(:)', sLen, 1);
        %primal residual r0
        curR0 = curZ0 - preY;
        curU0 = curU0 + curRho*curR0;
        %primal residual r1
        curR1 = curZ1 - curW;
        curU1 = curU1 + curRho*curR1;
        %primal residual r2
        diffW = cRblk*curW';
        curR2 = curZ2 - diffW;
        curU2 = curU2 + curRho*curR2;
        %compute the primal residual all
        rf = norm( curR0(:) ) + norm( curR1(:) ) + norm( curR2(:) );
        pDim = length( curR0(:) ) + length( curR1(:) ) + length( curR2(:) );
        
        %compute the dual variable to see if it is converge
        %dual residual s0
        curS0 = curRho*( D'*( curZ0 - preZ0 ) );
        curS1 = curRho*( curZ1 - preZ1 );
        curS2 = curRho*( cRblk'*( curZ2 - preZ2 ) );
        sf = norm( curS0(:) ) + norm( curS1(:) ) + norm( curS2(:) );
        nDim = length( curS0(:) ) + length( curS1(:) ) + length( curS2(:) );
        
        EPS_ABS = 1e-2;
        EPS_REL = 1e-4;
        
        [relEpsPri, idx] = max( [norm( curZ0(:) ) norm( preY(:) ) ...
            norm( curZ1(:) ), norm( curW(:) )...
            norm( curZ2(:) ), norm( diffW(:) )] );
        epsPri = sqrt(pDim)*EPS_ABS + relEpsPri*EPS_REL;
        
        tmp = D'*curU0;
        tmp = [tmp; sum(curU0)];
        if isempty( tmp )
            tmp = 0;
        end
        tmp2 = cRblk'*curU2;
        maxRel = max( [ norm( tmp(:) ), norm( u1(:) ), norm( tmp2(:) ) ] );
        epsDual = sqrt(nDim)*EPS_ABS + maxRel*EPS_REL;        
        resRecAry(itNumADMM, :, j) = [rf, epsPri, sf, epsDual];
        tmp = full( max( abs( curW(:) -  preW(:) ) ) );
        if isempty(tmp)
            tmp = 0;
        end
        if isempty(W)
            maxW = inf;
        else
            maxW = max(curW(:));
        end
%         time = toc;
%         fprintf( 'time: %g\n', time );        
        %% breaking condition
        if ( rf < epsPri && sf < epsDual ) && ( tmp < maxW*Wtol )
                %( max( abs(rf(:)) ) < 1e-3 && max( abs(sf(:)) ) < 1e-3 ) || ...
            break;
        end
        %% update the next stage rho
        if rf > 10*sf
            curRhoAry(itNumADMM+1, j) = 2 * curRhoAry(itNumADMM, j);
        elseif sf > 10*rf
            curRhoAry(itNumADMM+1, j) = 0.5 * curRhoAry(itNumADMM, j);
        else
            curRhoAry(itNumADMM+1, j) = curRhoAry(itNumADMM, j);
        end
        
        %% return the variables
        W(:, loc) = curW;
        W0(loc) = curW0;
        z0(:, loc) = curZ0;
        u0(:, loc) = curU0;
        z1(:, loc) = curZ1;
        u1(:, loc) = curU1;
        z2{j} = curZ2;
        u2{j} = curU2;
        LPAryADMM(itNumADMM+1, j) = LP_DL_Poiss( LINK_FUNC, aMatrix, Y, W, W0, D, lambda, 0, theta, scaleFactor, logFY, 0 );
        if j == BlkDS.blkNum
            if WHistFlag == 1
                WHistCell{itNumADMM} = W;
            end
        end
        fprintf( '(%d: %g) %g %g %g %g %g %g\n',itNumADMM, LPAryADMM(itNumADMM+1, j), rf, epsPri, sf, epsDual, curRhoAry(itNumADMM, j), full(tmp) );
        if itNumADMM == itNum
            if LPAryADMM(itNumADMM+1, j) >= LPAryADMM(1, j)
                if itNum <= oITNum*2
                    itNum = itNum + 5;
                    tmp = LPAryADMM;
                    LPAryADMM = zeros( itNum+1, BlkDS.blkNum );
                    LPAryADMM(1:(end-5),:) = tmp;
                else
                    stuckFlag = 1;
                    break;
                end
            end
        end
        itNumADMM = itNumADMM + 1;
    end
    if stuckFlag == 1
        break;
    end
end
initVar.W(rMIdx, :) = W; WResStruct.W = sparse( initVar.W );
initVar.W0 = W0; WResStruct.W0 = sparse( initVar.W0 );
initVar.z0 = z0; WResStruct.z0 = sparse( initVar.z0(:, :) );
initVar.z1(rMIdx, :) = z1; WResStruct.z1 = sparse( initVar.z1(:, :) );
for i = 1:BlkDS.blkNum
    initVar.z2{i}(:, rMIdx) = z2{i};
end
WResStruct.z2 = initVar.z2;
initVar.u0 = u0; WResStruct.u0 = sparse( initVar.u0(:, :) );
initVar.u1(rMIdx, :) = u1; WResStruct.u1 = sparse( initVar.u1(:, :) );
for i = 1:BlkDS.blkNum
    initVar.u2{i}(:, rMIdx) = u2{i};
end
WResStruct.u2 = initVar.u2;
WResStruct.LPAry = LPAryADMM;
WResStruct.rhoAry = curRhoAry;
WResStruct.resAry = resRecAry;
WResStruct.WDiff = tmp;
% WResStruct.kappa = kappa;
WResStruct.stuckFlag = stuckFlag;
if ~isempty(varargin) %WHistFlag == 1
    if varargin{1} == 1
        WResStruct.WHistCell = WHistCell;
    end
end
end

function C = replicateC( inD, repNum )
[i, j] = find( inD ~= 0 );
[sLen, mLen] = size(inD);
valVec = nonzeros( inD );
valVec = repmat( valVec, repNum, 1 );
yVec = zeros(length(i)*repNum, 1);
xVec = zeros(length(i)*repNum, 1);
unitLen = length(i);
for it = 1:repNum
    xVec( ((it-1)*unitLen+1):(unitLen*it) ) = j + ((it-1)*mLen);
    yVec( ((it-1)*unitLen+1):(unitLen*it) ) = i + ((it-1)*sLen);
end
C = sparse( yVec, xVec, valVec, sLen*repNum, mLen*repNum );

% for it = 2:repNum
%     tmp = S( ((it-1)*sLen+1):(it*sLen), : ) - S(1:sLen, :);
%     if ~isempty( find( tmp ~= 0, 1 ) )
%         fprintf( '%d\n', it );
%     end
% end
end

function C = rowCombSparseMatrix( A, B )
    [i, j] = find( A ~= 0 );
    [sLenA, mLenA] = size(A);
    [q, r] = find( B ~= 0 );
    [sLenB, mLenB] = size(B);
    if mLenA ~= mLenB
        fprintf( 'row combine sparse matrix error.\n' );
        C = [];
    else
        valVecA = nonzeros( A );
        valVecB = nonzeros( B );
        allVec = [valVecA; valVecB];
        q = q + sLenA;
        allYAxis = [i; q];
        allXAxis = [j; r];
        C = sparse( allYAxis, allXAxis, allVec, sLenA+sLenB, mLenA );
    end
end

function [C, target] = modifyCWithSparsity( CforUW, newWInfo, loc, mLen, nrPart1 )
    for i = 1:length(loc)
        cLoc = loc(i);
        cW = newWInfo{cLoc};
        cW = cW+(i-1)*mLen;
        cWRange = (1+(i-1)*mLen):(i*mLen);
        target = setdiff(cWRange, cW);
        CforUW(1:nrPart1, target) = 0;
    end
    [a, b] = find( CforUW((nrPart1+1):size(CforUW, 1), :)' ~= 0 );
    target = [];
    for i = 1:2:length(a)
        [q, r] = ind2sub([mLen, length(loc)], a(i:(i+1)) );
        if isempty( find( newWInfo{loc(r(1))} == q(1), 1 ) ) || isempty( find( newWInfo{loc(r(2))} == q(2), 1 ) )
            target = [ target b(i)];
        end
    end
    target = target + nrPart1;
    CforUW( target, :) = 0;
    C = CforUW;
end
