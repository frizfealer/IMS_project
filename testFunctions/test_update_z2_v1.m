function [ output_args ] = test_update_z2_v1( testCaseNum )
%test function for update_z1_v2
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
    rho = 1;
    W = simData.gW(:, :);
    BlkDS = conBLKDS( simData.gY );
    [Rall] = genSparseGroupingMatrix( HEIGHT, WIDTH, 1 );
    for i = 1:BlkDS.blkNum
        Rblk{i} = [];
        Rb = Rall(:, BlkDS.B2GMap{i});
        ins = sum(abs(Rb), 2);
        Rb(ins<2,:) = [];
        Rblk{i} = Rb;
    end
    fVal = zeros( 2, BlkDS.blkNum );
    for i = 1:BlkDS.blkNum
        loc = BlkDS.B2GMap{i};
        curW = W(:, loc);
        curR = Rblk{i};
        curZ2 = zeros( size(curR, 1), MLEN );
        preZ2 = curZ2;
        curU2 = zeros( size(curZ2) );
        L1Flag = 1;
        j = 1e-1;
        RMat = genSparseGroupingMatrix2( curR, MLEN, 0 );
        fVal(1, i) = Z2_term_func( curR, curW, preZ2, curU2, rho, j, L1Flag );
        [ curZ2 ] = updatez2_v1( curZ2, curU2, rho, curR, curW, j, L1Flag );
        fVal(2, i) = Z2_term_func( curR, curW, curZ2, curU2, rho, j, L1Flag );
        gZ2 = curR*curW';
%         figure; imagesc( abs( curZ2-gZ2 ) ); colorbar;
    end
    fVal = zeros( 4, 100 );
    dfVal = zeros( 4, 100 );
    theta = 1e-1;
    for j = 1:100
        for i = 1:BlkDS.blkNum
            loc = BlkDS.B2GMap{i};
            curW = W(:, loc);
            curR = Rblk{i};
            curZ2 = curR*curW';
            preZ2 = curZ2;
            curU2 = zeros( size(curZ2) );
            L1Flag = 1;
            fVal((i-1)*2+1, j) = Z2_term_func( curR, curW, preZ2, curU2, rho, theta, L1Flag );
            dfVal((i-1)*2+1, j) = length( find( preZ2(:) >= 1e-3 ) );
            [ curZ2 ] = updatez2_v1( curZ2, curU2, rho, curR, curW, theta, L1Flag );
            fVal((i-1)*2+2, j) = Z2_term_func( curR, curW, curZ2, curU2, rho, theta, L1Flag );
            dfVal((i-1)*2+2, j) = length( find( curZ2(:) >= 1e-3 ) );
        end
        theta = 1e-1 * (j+1);
        if j == 5 && L1Flag == 1
            figure; hist( curZ2(:), 100 );
        end
    end
   figure; subplot(1, 2, 1); plot(dfVal(1:2, :)' ); subplot(1, 2, 2); plot(dfVal(3:4, :)' );
   figure; subplot(1, 2, 1); plot(fVal(1:2, :)' ); subplot(1, 2, 2); plot(fVal(3:4, :)' );
elseif testCaseNum == 2
    load( 'D:\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LN\100831_348_136_12_40hr_0_1XLB_LN_input.mat' );
    load( 'D:\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LN\exp_100831_348_136_12_40hr_0_1XLB_LN_4_res.mat' );
    [sLen, hei, wid] = size(dataCube);
    W = expRec.outW;
    fVal = zeros(4, 1);
    dfVal = zeros( 2, 1 );
    [Rall] = genSparseGroupingMatrix( hei, wid, 1 );
    for i = 1:BlkDS.blkNum
        Rblk{i} = [];
        Rb = Rall(:, BlkDS.B2GMap{i});
        ins = sum(abs(Rb), 2);
        Rb(ins<2,:) = [];
        Rblk{i} = Rb;
    end
    theta = 1;
    for j = 1:BlkDS.blkNum
        loc = BlkDS.B2GMap{j};
        curW = W(:, loc);
        curR = Rblk{j};
        curZ2 = curR*curW';
        preZ2 = curZ2;
        curU2 = zeros( size( curZ2 ) );
        L1Flag = 1;
        rho = 1;
        fVal((j-1)*2+1, 1) = Z2_term_func( curR, curW, preZ2, curU2, rho, theta, L1Flag );
        dfVal((j-1)*2+1, 1) = length( find( preZ2(:) >= 1e-3 ) );
        [ curZ2 ] = updatez2_v1( curZ2, curU2, rho, curR, curW, theta, L1Flag );
        fVal((j-1)*2+2, 1) = Z2_term_func( curR, curW, curZ2, curU2, rho, theta, L1Flag );
        dfVal((j-1)*2+2, 1) = length( find( curZ2(:) >= 1e-3 ) );
    end
    fVal
    dfVal
    fVal = zeros(4, 100 );
    dfVal = zeros( 4, 100 );
    theta = 1e-2;
    for i = 1:100
        fprintf( '%d ', i );
        for j = 1:BlkDS.blkNum
            loc = BlkDS.B2GMap{j};
            curW = W(:, loc);
            curR = Rblk{j};
            curZ2 = curR*curW';
            preZ2 = curZ2;
            curU2 = zeros( size( curZ2 ) );
            L1Flag = 0;
            rho = 1;
            fVal((j-1)*2+1, i) = Z2_term_func( curR, curW, preZ2, curU2, rho, theta, L1Flag );
            dfVal((j-1)*2+1, i) = length( find( preZ2(:) >= 1e-3 ) );
            [ curZ2 ] = updatez2_v1( curZ2, curU2, rho, curR, curW, theta, L1Flag );
            fVal((j-1)*2+2, i) = Z2_term_func( curR, curW, curZ2, curU2, rho, theta, L1Flag );
            dfVal((j-1)*2+2, i) = length( find( curZ2(:) >= 1e-3 ) );
        end
        theta = 1e-2*(i+1);
    end
    figure; subplot(1, 2, 1); plot(dfVal(1:2, :)' ); subplot(1, 2, 2); plot(dfVal(3:4, :)' );
    figure; subplot(1, 2, 1); plot(fVal(1:2, :)' ); subplot(1, 2, 2); plot(fVal(3:4, :)' );
end

end

function [ val ] = Z2_term_func( R, W, Z2, U2, rho, theta, L1Flag )
%this function should be minimized.
val1 = R*W'-1/rho*U2-Z2;
val1 = (norm(val1(:))^2)*rho/2;
if L1Flag == 1
    val2 = norm( Z2(:), 1 )*theta;
else
    val2 = (norm( Z2(:) )^2)*theta;
end
val = val1+val2;
end

