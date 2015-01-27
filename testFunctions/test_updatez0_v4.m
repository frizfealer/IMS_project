function [  ] = test_updatez0_v4( testCaseNum )
%test_updatez0_v4 test function for updatez0_v4
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
    LINK_FUNC = 'log';
    [simData] = synthesizeData_Poisson( LINK_FUNC, SLEN, MLEN, HEIGHT, WIDTH, [], DOptions, WOptions, verbose );
    %% test on updatez0_v4 function, no held out data
    alpha = ones(HEIGHT, WIDTH);
    inY = simData.gY;
    BlkDS = conBLKDS(inY);
    gD = simData.gD;
    gW = simData.gW;
    gW0 = simData.gW0;
    U0 = zeros( SLEN, HEIGHT*WIDTH );
    curRho = 1;
    scaleFactor = 1e-8;
    Z0 = zeros( SLEN, HEIGHT*WIDTH );
    funcVal(1) = Z0_termFunc( inY, Z0, curRho, gD, gW, gW0, U0, scaleFactor );
    for j = 1:BlkDS.blkNum
            loc = BlkDS.B2GMap{j};
            locAlpha1 = intersect( loc, find( alpha(:) == 1 ) );
            %update for training data
            curY = inY(:, locAlpha1);
            curW = gW(:, locAlpha1);
            curW0 = gW0(locAlpha1);
            curZ0 = Z0(:, locAlpha1);
            curU0 = U0(:, locAlpha1);
            [ tmpZ0 ] = updatez0_v4( 1, curY, curZ0, gD, curW, curW0, curU0, curRho, scaleFactor );
            Z0(:, locAlpha1) = tmpZ0;
            funcVal(j+1) = Z0_termFunc( inY, Z0, curRho, gD, gW, gW0, U0, scaleFactor );
    end
    funcVal
    statMap = zeros( SLEN, HEIGHT*WIDTH );
    for j = 1:BlkDS.blkNum
        loc = BlkDS.B2GMap{j};
        tmp = abs( (Z0(:, loc)) - log(inY(:, loc)) );
        statMap(:, loc) = tmp;
    end
    max(max(statMap))
    figure;
    for i = 1:SLEN
        subplot( 2, 5, i ); imagesc( reshape( statMap(i, :), HEIGHT, WIDTH ) );
    end
elseif testCaseNum == 2
    load( 'D:\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LN\100831_348_136_12_40hr_0_1XLB_LN_input.mat' );
    load( 'D:\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LN\exp_100831_348_136_12_40hr_0_1XLB_LN_4_res.mat' );
    method = {'csd', 'bb', 'pnewton0', 'scg', 'pcg', 'lbfgs', 'newton'};
    for i = 1:7
        fprintf( '%s\n', method{i});
        inY = dataCube(:, :);
        Z0 = zeros( size( expRec.z0 ) );
        curRho = 1;
        gD = expRec.outD;
        gW = expRec.outW;
        gW0 = expRec.outW0;
        U0 = zeros( size(inY, 1), size(inY, 2)*size(inY, 3) );
        scaleFactor = 1;
        funcVal(1) = Z0_termFunc( inY, Z0, curRho, gD, gW, gW0, U0, scaleFactor );
        alpha = ones( size(inY, 2), size(inY, 3) );
        tic
        for j = 1:BlkDS.blkNum
            loc = BlkDS.B2GMap{j};
            locAlpha1 = intersect( loc, find( alpha(:) == 1 ) );
            %update for training data
            curY = inY(:, locAlpha1);
            curW = gW(:, locAlpha1);
            curW0 = gW0(locAlpha1);
            curZ0 = Z0(:, locAlpha1);
            curU0 = U0(:, locAlpha1);
            [ tmpZ0 ] = updatez0_v4( 1, curY, curZ0, gD, curW, curW0, curU0, curRho, scaleFactor, 1000, method{i} );
            Z0(:, locAlpha1) = tmpZ0;
            funcVal(j+1) = Z0_termFunc( inY, Z0, curRho, gD, gW, gW0, U0, scaleFactor );
        end
        toc
        funcVal
    end
elseif testCaseNum == 3
    load( 'D:\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LN\100831_348_136_12_40hr_0_1XLB_LN_input.mat' );
    load( 'D:\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LN\exp_100831_348_136_12_40hr_0_1XLB_LN_4_res.mat' );
    inY = dataCube(:, :);
    Z0 = zeros( size( expRec.z0 ) );
    curRho = 1;
    gD = expRec.outD;
    gW = expRec.outW;
    gW0 = expRec.outW0;
    U0 = zeros( size(inY, 1), size(inY, 2)*size(inY, 3) );
    scaleFactor = 1;
    funcVal(1) = Z0_termFunc( inY, Z0, curRho, gD, gW, gW0, U0, scaleFactor );
    alpha = ones( size(inY, 2), size(inY, 3) );
    tic
    for i = 1:size(inY, 2)*size(inY, 3)
        [ Z0(:, i) ] = updatez0_v3( alpha(i), inY(:,i), Z0(:,i), gD, gW(:, i), gW0(i), U0(:,i), curRho, scaleFactor );
    end
    toc
    funcVal(2) = Z0_termFunc( inY, Z0, curRho, gD, gW, gW0, U0, scaleFactor );
elseif testCaseNum == 4
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
    LINK_FUNC = 'identity';
    [simData] = synthesizeData_Poisson( LINK_FUNC, SLEN, MLEN, HEIGHT, WIDTH, [], DOptions, WOptions, verbose );
    %% test on updatez0_v4 function, no held out data
    alpha = ones(HEIGHT, WIDTH);
    inY = simData.gY;
    BlkDS = conBLKDS(inY);
    gD = simData.gD;
    gW = simData.gW;
    gW0 = simData.gW0;
    U0 = zeros( SLEN, HEIGHT*WIDTH );
    curRho = 1;
    scaleFactor = 1;
    Z0 = abs(randn( SLEN, HEIGHT*WIDTH )) + 1e-32;
    funcVal(1) = Z0_termFunc_identity( inY, Z0, curRho, gD, gW, gW0, U0, scaleFactor );
    for j = 1:BlkDS.blkNum
        loc = BlkDS.B2GMap{j};
        locAlpha1 = intersect( loc, find( alpha(:) == 1 ) );
        %update for training data
        curY = inY(:, locAlpha1);
        curW = gW(:, locAlpha1);
        curW0 = gW0(locAlpha1);
        curZ0 = Z0(:, locAlpha1);
        curU0 = U0(:, locAlpha1);
        [ tmpZ0 ] = updatez0_v4( 1, LINK_FUNC, curY, curZ0, gD, curW, curW0, curU0, curRho, scaleFactor, 200, 'newton' );
        Z0(:, locAlpha1) = tmpZ0;
        funcVal(j+1) = Z0_termFunc_identity( inY, Z0, curRho, gD, gW, gW0, U0, scaleFactor );
    end
    funcVal
    statMap = zeros( SLEN, HEIGHT*WIDTH );
    for j = 1:BlkDS.blkNum
        loc = BlkDS.B2GMap{j};
        tmp = abs( (Z0(:, loc)) - inY(:, loc) );
        statMap(:, loc) = tmp;
    end
    max(max(statMap))
    figure;
    for i = 1:SLEN
        subplot( 2, 5, i ); imagesc( reshape( statMap(i, :), HEIGHT, WIDTH ) );
    end
elseif testCaseNum == 5
    %% synthesize data first
    SLEN = 10;
    MLEN = 10;
    HEIGHT = 30;
    WIDTH = 30;
    DOptions.sparsePrec = 1;
    DOptions.coheMax = 1;
    WOptions.supPrec = 1;
    WOptions.sparsePrec = 1;
    verbose = 0;
    LINK_FUNC = 'log_gaussain';
    [simData] = synthesizeData_Poisson( LINK_FUNC, 'L1', SLEN, MLEN, HEIGHT, WIDTH, [], DOptions, WOptions, verbose );
    %% test on updatez0_v4 function, no held out data
    alpha = ones(HEIGHT, WIDTH);
    inY = simData.gY;
    BlkDS = conBLKDS(inY);
    gD = simData.gD;
    gW = simData.gW;
    gW0 = simData.gW0;
    U0 = zeros( SLEN, HEIGHT*WIDTH );
    curRho = 1;
    scaleFactor = 1;
%     Z0 = abs(randn( SLEN, HEIGHT*WIDTH )) + 1e-32;
    Z0 = simData.gD*gW(:,:);
    funcVal(1) = Z0_termFunc_log_gaussain( inY, Z0, curRho, gD, gW, gW0, U0, scaleFactor );
    for j = 1:BlkDS.blkNum
        loc = BlkDS.B2GMap{j};
        locAlpha1 = intersect( loc, find( alpha(:) == 1 ) );
        %update for training data
        curY = inY(:, locAlpha1);
        curW = gW(:, locAlpha1);
        curW0 = gW0(locAlpha1);
        curZ0 = Z0(:, locAlpha1);
        curU0 = U0(:, locAlpha1);
        [ tmpZ0 ] = updatez0_v4( 1, LINK_FUNC, curY, curZ0, gD, curW, curW0, curU0, curRho, scaleFactor, 200, 'newton' );
        Z0(:, locAlpha1) = tmpZ0;
        funcVal(j+1) = Z0_termFunc_log_gaussain( inY, Z0, curRho, gD, gW, gW0, U0, scaleFactor );
    end
    funcVal
    statMap = zeros( SLEN, HEIGHT*WIDTH );
    preY = simData.gD*gW(:, :);
    for j = 1:BlkDS.blkNum
        loc = BlkDS.B2GMap{j};
        tmp = abs( (Z0(:, loc)) - preY(:, loc) );
        statMap(:, loc) = tmp;
    end
    max(max(statMap))
    figure;
    for i = 1:SLEN
        subplot( 2, 5, i ); imagesc( reshape( statMap(i, :), HEIGHT, WIDTH ) );
    end
end
end
function [ val ] = Z0_termFunc( inY, Z0, rho, gD, gW, gW0, U0, scaleFac )
res1 = -gD*gW(:, :) - repmat( gW0(:)', size(gD, 1), 1 ) + 1/rho*U0;
eZ0 = exp(Z0);
%     fprintf( 'll = %g, pp = %g\n', full(scaleFac*sum( ( -y.*z0 + eZ0 ))), full(rho / 2 * sum( (z0 +res1).^2 )));
val = sum( sum( scaleFac*( -inY(:, :).*Z0 + eZ0 ) ) ) + rho / 2 * sum( sum( (Z0 +res1).^2 ) );
end

function [ val ] = Z0_termFunc_identity( inY, Z0, rho, gD, gW, gW0, U0, scaleFac )
res1 = -gD*gW(:, :) - repmat( gW0(:)', size(gD, 1), 1 ) + 1/rho*U0;
%     fprintf( 'll = %g, pp = %g\n', full(scaleFac*sum( ( -y.*z0 + eZ0 ))), full(rho / 2 * sum( (z0 +res1).^2 )));
val = sum( sum( scaleFac*( -inY(:, :).*log(Z0+1e-32) + Z0 ) ) ) + rho / 2 * sum( sum( (Z0 +res1).^2 ) );
end

function [ val ] = Z0_termFunc_log_gaussain( inY, Z0, rho, gD, gW, gW0, U0, scaleFac )
res1 = -gD*gW(:, :) - repmat( gW0(:)', size(gD, 1), 1 ) + 1/rho*U0;
%     fprintf( 'll = %g, pp = %g\n', full(scaleFac*sum( ( -y.*z0 + eZ0 ))), full(rho / 2 * sum( (z0 +res1).^2 )));
val = sum( sum( scaleFac*(log(inY(:, :)+1e-32) - Z0).^2 ) ) + rho / 2 * sum( sum( (Z0 +res1).^2 ) );
end