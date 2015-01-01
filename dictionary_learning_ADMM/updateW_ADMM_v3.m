function [WResStruct] = updateW_ADMM_v3( Y, D, aMatrix, itNum, lambda, theta, L1Flag, logFY, initVar, varargin )

%% construct metaData
BlkDS = conBLKDS( Y );
[sLen, hei, wid] = size( Y );
[Rall] = genSparseGroupingMatrix( hei, wid, 1 );
for i = 1:BlkDS.blkNum
    Rblk{i} = [];
    Rb = Rall(:, BlkDS.B2GMap{i});
    ins = sum(abs(Rb), 2);
    Rb(ins<2,:) = [];
    Rblk{i} = Rb;
end

%% variable setting
fprintf( 'iteration number for W_update_ADMM = %d\n', itNum );
temp = (aMatrix==0)&(BlkDS.indMap==1);
fprintf( 'held-out sample number = %d, training sample number = %d\n', length( find( temp == 1 ) ), length( find( aMatrix == 1 ) ) );
fprintf( 'lambda = %d, theta = %d\n', lambda, theta );
mLen = size(D, 2);
%% initializing z0 as log(inY)

if isempty( initVar )
    z0 = log(Y); z0(z0==-inf)=0; %initialize zo as log(Y)
    %if not in traing set, we shoud not initialize z0 according to it.
    %take the average values
    idx = intersect( find(aMatrix==1), find(BlkDS.indMap==1) );
    tmp = sum( z0(:, idx), 2 ) / ( length( idx ) );
    for i = 1:hei*wid
        %if is in data area, but not in training set,
        %set it to average spectrum
        if BlkDS.indMap(i) && ~aMatrix(i)
            z0(:, i) = tmp;
        end
    end
    z1 = zeros( mLen, wid*hei ); 
    z2 = cell( BlkDS.blkNum, 1 );
    for i = 1:BlkDS.blkNum
        z2{i} = zeros( length(Rblk{i}), mLen );
    end
    W = zeros( mLen, wid*hei ); W0 = zeros( hei, wid );
else
    z0 = initVar.z0(:, :); z1 = initVar.z1(:, :); z2 = initVar.z2;
    W = initVar.W(:, :); W0 = initVar.W0(:, :);
end
scaleFactor =  1 / ( max( Y(:) ) / max( z0(:) ) ) * 100;

LPAryADMM = zeros( itNum+1, BlkDS.blkNum );
resRecAry = zeros( itNum, 4 );
curRhoAry = zeros( itNum, 1 );
curRhoAry(1) = 1;
u0 = sparse( zeros( size( z0(:, :) ) ) );
u1 = sparse( zeros( size( z1(:, :) ) ) );
u2 = cell( BlkDS.blkNum, 1 );
for i = 1:BlkDS.blkNum
    u2{i} = sparse( zeros( size( z2{i}(:, :) ) ) );
end

moutD = [D ones(sLen, 1)];
tmpEye = eye(mLen+1, mLen+1);
tmpEye(end, end) = 0;
aPartofCforUW = sparse( [moutD; tmpEye] );
%% main optimization process
for j = 1:BlkDS.blkNum
    LPAryADMM(1, j) = LP_DL_Poiss( aMatrix, Y, W, W0, D, lambda, 0, theta, scaleFactor, logFY, 0 );
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
    RforUW = genSparseGroupingMatrix2( cRblk, mLen, 1 );
    CforUW = rowCombSparseMatrix( CforUW, RforUW );
    for itNumADMM = 1:itNum
        tic;
        preW = curW;
        fprintf( 'LP: %g ', LPAryADMM(itNumADMM) );
        curRho = curRhoAry(itNumADMM);
        %% update W
        fprintf( 'updating W... ' );
        [ w_w0 ] = update_w_w0_v4( curZ0, curU0, curZ1, curU1, curZ2, curU2, CforUW, curRho, [curW; curW0(:)'] );
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
        [ curZ0(:, locAlpha0) ] = updatez0_v4( 0, curY(:, locAlpha0), curZ0(:, locAlpha0), D, ...
            curW(:, locAlpha0), curW0(locAlpha0), curU0(:, locAlpha0), curRho, scaleFactor );
        %update for training data
        [ curZ0(:, locAlpha1) ] = updatez0_v4( 1, curY(:, locAlpha1), curZ0(:, locAlpha1), D, ...
            curW(:, locAlpha1), curW0(locAlpha1), curU0(:, locAlpha1), curRho, scaleFactor );
        %% update z1
        fprintf( 'updating z1... ' );
        preZ1 = curZ1;
        [ curZ1 ] = updatez1_v2( curZ1, curU1, curW, curRho, lambda );
        %% update z2
        fprintf( 'updating z2... ' );
        preZ2 = curZ2;
        [ curZ2 ] = updatez2_v1( curZ2, curU2, curRho, cRblk, curW, theta, L1Flag );
        %% update u0, u1
        fprintf( 'updating u0, u1, u2... ' );
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
        
        EPS_ABS = 1e-4;
        EPS_REL = 1e-3;
        
        relEpsPri = max( [norm( curZ0(:) ) norm( preY(:) ) ...
            norm( curZ1(:) ), norm( curW(:) )...
            norm( curZ2(:) ), norm( diffW(:) )] );
        epsPri = sqrt(pDim)*EPS_ABS + relEpsPri*EPS_REL;
        
        tmp = D'*curU0;
        tmp = [tmp; sum(curU0)];
        tmp2 = cRblk'*curU2;
        maxRel = max( [ norm( tmp(:) ), norm( u1(:) ), norm( tmp2(:) ) ] );
        epsDual = sqrt(nDim)*EPS_ABS + maxRel*EPS_REL;
        
        resRecAry(itNumADMM, :) = [rf, epsPri, sf, epsDual];
        tmp = full( max( abs( curW(:) -  preW(:) ) ) );
        fprintf( '%d: %g %g %g %g %g %g ',itNumADMM, rf, epsPri, sf, epsDual, curRhoAry(itNumADMM), full(tmp) );
        time = toc;
        fprintf( 'time: %g\n', time );        
        %% breaking condition
        if ( rf < epsPri && sf < epsDual ) || ( tmp < 1e-3 )
                %( max( abs(rf(:)) ) < 1e-3 && max( abs(sf(:)) ) < 1e-3 ) || ...
            break;
        end
        %% update the next stage rho
        if rf > 10*sf
            curRhoAry(itNumADMM+1) = 2 * curRhoAry(itNumADMM);
        elseif sf > 10*rf
            curRhoAry(itNumADMM+1) = 0.5 * curRhoAry(itNumADMM);
        else
            curRhoAry(itNumADMM+1) = curRhoAry(itNumADMM);
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
        LPAryADMM(itNumADMM+1, j) = LP_DL_Poiss( aMatrix, Y, W, W0, D, lambda, 0, theta, scaleFactor, logFY, 0 );
        
        if j == BlkDS.blkNum
            if ~isempty(varargin) %WHistFlag == 1
                if varargin{1} == 1
                    WHistCell{itNumADMM} = W;
                end
            end
        end

    end
end
WResStruct.W = sparse( W );
WResStruct.W0 = sparse( W0 );
WResStruct.z0 = sparse( z0(:, :) );
WResStruct.z1 = sparse( z1(:, :) );
WResStruct.z2 = z2;
WResStruct.u0 = sparse(u0(:, :) );
WResStruct.u1 = sparse(u1(:, :) );
WResStruct.u2 = u2;
WResStruct.LPAry = sparse( LPAryADMM );
WResStruct.rhoAry = sparse( curRhoAry );
WResStruct.resAry = sparse( resRecAry );
WResStruct.WDiff = sparse( tmp );
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
