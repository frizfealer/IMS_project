function [  ] = test_update_w_w0_v4( testCaseNum )
%test_updatez0_v4 test function for update_w_w0_v4
if testCaseNum == 1
    %% synthesize data first
    SLEN = 10;
    MLEN = 10;
    HEIGHT = 20;
    WIDTH = 20;
    DOptions.sparsePrec = 1;
    DOptions.coheMax = 1;
    WOptions.supPrec = 1;
    WOptions.sparsePrec = 1;
    verbose = 0;
    [simData] = synthesizeData_Poisson( SLEN, MLEN, HEIGHT, WIDTH, [], DOptions, WOptions, verbose );
    %% test on update_w_w0_v4 function, no held out data
    alpha = ones(HEIGHT, WIDTH);
    inY = simData.gY;
    BlkDS = conBLKDS(inY);
    gD = simData.gD;
    W = zeros( size( simData.gW ) );
    W0 = zeros( size( simData.gW0 ) );
    curRho = 1;
    Z0 = simData.gD*simData.gW(:, :);
    U0 = zeros( size( Z0 ) );
    Z1 = simData.gW;
    U1 = zeros( size( Z1 ) );
    [Rall] = genSparseGroupingMatrix( HEIGHT, WIDTH, 1 );
    for i = 1:BlkDS.blkNum
        Rblk{i} = [];
        Rb = Rall(:, BlkDS.B2GMap{i});
        ins = sum(abs(Rb), 2);
        Rb(ins<2,:) = [];
        Rblk{i} = Rb;
    end
    moutD = [gD ones(SLEN, 1)];
    tmpEye = eye(MLEN+1, MLEN+1);
    tmpEye(end, end) = 0;
    partofCforUW = sparse( [moutD; tmpEye] );
    
    for j = 1:BlkDS.blkNum
            loc = BlkDS.B2GMap{j};
            %update for training data
            curW = W(:, loc);
            curW0 = W0(loc);
            curZ0 = Z0(:, loc);
            curU0 = U0(:, loc);
            curZ1 = Z1(:, loc);
            curU1 = U1(:, loc);
            curZ2 = Rblk{j}*simData.gW(:,loc)';
            curU2 = zeros( size(curZ2) );
            funcVal(1) = w_w0_term_func( curZ0, curZ1, curZ2, curRho, gD, Rblk{j}, curW, curW0, curU0, curU1, curU2 );
            CforUW = replicateC( partofCforUW, length(loc) );
            RforUW = genSparseGroupingMatrix2( Rblk{j}, MLEN, 1 );
            CforUW = rowCombSparseMatrix( CforUW, RforUW );
            [ w_w0 ] = update_w_w0_v4( curZ0, curU0, curZ1, curU1, curZ2, curU2, CforUW, curRho, [curW; curW0(:)'] );
            tmp = reshape(w_w0, MLEN+1, length(loc));
            curW = tmp(1:(end-1), :);
            curW0 = tmp(end,:);
            funcVal(2) = w_w0_term_func( curZ0, curZ1, curZ2, curRho, gD, Rblk{j}, curW, curW0, curU0, curU1, curU2 );
            W(:, loc) = curW;
            W0(loc) = curW0;
    end
    figure; 
    for i = 1:MLEN
        tmp = abs(simData.gD*W(:,:)+repmat(W0(:)', SLEN, 1)-simData.gD*simData.gW(:,:));
        subplot(2, 5, i); imagesc( reshape( tmp(i,:),20, 20) ); colorbar;
    end
    figure;
    for i = 1:5
        subplot(2, 5, i); imagesc( reshape(simData.gW(i,:), 20, 20) ); colorbar;
        subplot(2, 5, i+5); imagesc( reshape(W(i,:), 20, 20) ); colorbar;
    end
    figure;
    for i = 1:5
        subplot(2, 5, i); imagesc( reshape(simData.gW(i+5,:), 20, 20) ); colorbar;
        subplot(2, 5, i+5); imagesc( reshape(W(i+5,:), 20, 20) ); colorbar;
    end
elseif testCaseNum == 2
    load( 'D:\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LN\100831_348_136_12_40hr_0_1XLB_LN_input.mat' );
    load( 'D:\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LN\exp_100831_348_136_12_40hr_0_1XLB_LN_4_res.mat' );
        inY = dataCube(:, :);
        [sLen, hei, wid] = size(dataCube);
        Z0 = expRec.z0;
        U0 = zeros( size( Z0 ) );
        Z1 = expRec.z1;
        U1 = zeros( size( Z1 ) );
        [Rall] = genSparseGroupingMatrix( hei, wid, 1 );
        for i = 1:BlkDS.blkNum
            Rblk{i} = [];
            Rb = Rall(:, BlkDS.B2GMap{i});
            ins = sum(abs(Rb), 2);
            Rb(ins<2,:) = [];
            Rblk{i} = Rb;
        end
        curRho = 1;
        gD = expRec.outD;
        [~, mLen] = size(gD);
        W = expRec.outW;
        W0 = expRec.outW0;
        scaleFactor = 1;
        moutD = [gD ones(sLen, 1)];
        tmpEye = eye(mLen+1, mLen+1);
        tmpEye(end, end) = 0;
        partofCforUW = sparse( [moutD; tmpEye] );
        for j = 1:BlkDS.blkNum
            loc = BlkDS.B2GMap{j};
            %update for training data
            curW = W(:, loc);
            curW0 = W0(loc);
            curZ0 = Z0(:, loc);
            curU0 = U0(:, loc);
            curZ1 = Z1(:, loc);
            curU1 = U1(:, loc);
            curZ2 = Rblk{j}*curW';
            curU2 = zeros( size(curZ2) );
            funcVal(1) = w_w0_term_func( curZ0, curZ1, curZ2, curRho, gD, Rblk{j}, curW, curW0, curU0, curU1, curU2 );
            CforUW = replicateC( partofCforUW, length(loc) );
            RforUW = genSparseGroupingMatrix2( Rblk{j}, mLen, 1 );
            CforUW = rowCombSparseMatrix( CforUW, RforUW );
            [ w_w0 ] = update_w_w0_v4( curZ0, curU0, curZ1, curU1, curZ2, curU2, CforUW, curRho, [curW; curW0(:)'] );
            tmp = reshape(w_w0, mLen+1, length(loc));
            curW = tmp(1:(end-1), :);
            curW0 = tmp(end, :);
            funcVal(2) = w_w0_term_func( curZ0, curZ1, curZ2, curRho, gD, Rblk{j}, curW, curW0, curU0, curU1, curU2 );
            W(:, loc) = curW;
            W0(loc) = curW0;
            funcVal
        end
end

end

function [ val ] = w_w0_term_func( Z0, Z1, Z2, rho, gD, R, W, W0, U0, U1, U2 )
%this function should be minimized.
val1 = ( Z0 - gD*W - repmat(W0(:)', size(Z0, 1), 1) + 1/rho*U0 );
val1 = norm( val1(:) );
val2 = ( Z1 - W + 1/rho*U1 );
val2 = norm( val2(:) );
val3 = ( Z2 - R*W' + 1/rho*U2 );
val3 = norm( val3(:) );
val = val1+val2+val3;
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