function [ output_args ] = test_updateW_ADMM_v3( testCaseNum )
%test_updateW_ADMM_v3 test updateW_ADMM_v3
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
    aMatrix = ones( HEIGHT, WIDTH );
    itNum = 100;
    lambda = 1e-2;
    theta = 1*1e-32;
    L1Flag = 0;
    logFY = [];
    initVar = [];
    [WResStruct] = updateW_ADMM_v3( simData.gY, simData.gD, aMatrix, itNum, lambda, theta, L1Flag, logFY, initVar );
    for j = 1:2
        figure;
        for i = 1:5
            cW = WResStruct.W;
            gW = simData.gW;
            subplot(2, 5, i); imagesc( reshape(cW(i*j, :), HEIGHT, WIDTH ) ); colorbar;
            subplot(2, 5, i+5); imagesc( reshape( gW(i*j, :), HEIGHT, WIDTH) ); colorbar;
        end
    end
    itNum = 200;
    [WResStruct] = updateW_ADMM_v3( simData.gY, simData.gD, aMatrix, itNum, lambda, theta, L1Flag, logFY, initVar );
    for j = 1:2
        figure;
        for i = 1:5
            cW = WResStruct.W;
            gW = simData.gW;
            subplot(2, 5, i); imagesc( reshape(cW(i*j, :), HEIGHT, WIDTH ) ); colorbar;
            subplot(2, 5, i+5); imagesc( reshape( gW(i*j, :), HEIGHT, WIDTH) ); colorbar;
        end
    end

elseif testCaseNum == 2
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
    itNum = 200;
    lambda = 1e-32;
    theta = 1e-2;
    L1Flag = 1;
    logFY = [];
    initVar = [];
    BlkDS = conBLKDS( simData.gY );
    [ aMatrix ] = leaveoutCBPatternsData( BlkDS, 10 );
    [ scaleFactor ] = computeScaleFactor( simData.gY, aMatrix );
    figure; imagesc(aMatrix);
    [WResStruct] = updateW_ADMM_v3( simData.gY, simData.gD, aMatrix, itNum, lambda, theta, L1Flag, logFY, initVar, scaleFactor );
    for j = 1:2
        figure;
        for i = 1:5
            cW = WResStruct.W;
            gW = simData.gW;
            subplot(2, 5, i); imagesc( reshape(cW(i*j, :), HEIGHT, WIDTH ) ); colorbar;
            subplot(2, 5, i+5); imagesc( reshape( gW(i*j, :), HEIGHT, WIDTH) ); colorbar;
        end
    end
    figure; subplot(1, 2, 1); plot(WResStruct.LPAry(:, 1)); subplot(1, 2, 2); plot(WResStruct.LPAry(:, 2) );
    
    aMatrix = BlkDS.indMap;
    figure; imagesc(aMatrix);
    [WResStruct] = updateW_ADMM_v3( simData.gY, simData.gD, aMatrix, itNum, lambda, theta, L1Flag, logFY, initVar, scaleFactor );
    for j = 1:2
        figure;
        for i = 1:5
            cW = WResStruct.W;
            gW = simData.gW;
            subplot(2, 5, i); imagesc( reshape(cW(i*j, :), HEIGHT, WIDTH ) ); colorbar;
            subplot(2, 5, i+5); imagesc( reshape( gW(i*j, :), HEIGHT, WIDTH) ); colorbar;
        end
    end
    figure; subplot(1, 2, 1); plot(WResStruct.LPAry(:, 1)); subplot(1, 2, 2); plot(WResStruct.LPAry(:, 2) );

elseif testCaseNum == 3
    load( 'D:\Google ¶³ºİµwºĞ\realExp_100831_348_136_12_40hr_0-1XLB_LN\100831_348_136_12_40hr_0_1XLB_LN_input.mat' );
    load( 'D:\Google ¶³ºİµwºĞ\realExp_100831_348_136_12_40hr_0-1XLB_LN\exp_100831_348_136_12_40hr_0_1XLB_LN_4_res.mat' );
    [~, hei, wid] = size(dataCube);
    aMatrix = ones( hei, wid );
    itNum = 100;
    lambda = 0;
    theta = 1e-32;
    L1Flag = 1;
    logFY = [];
    initVar = [];
    [ scaleFactor ] = computeScaleFactor( dataCube, aMatrix );
    iD = initD( dataCube, nDTemplate, 'NNMF', nDIonName );
    [WResStruct] = updateW_ADMM_v3( dataCube, iD, aMatrix, itNum, lambda, theta, L1Flag, logFY, initVar, scaleFactor, 0, 5*1e-2 );
    figure; subplot(1, 2, 1); plot(WResStruct.LPAry(:, 1)); subplot(1, 2, 2); plot(WResStruct.LPAry(:, 2) );
    ins = WResStruct.W(:, BlkDS.indMap == 1 );
    tmp = zeros( size(ins, 1), 1);
    for i = 1:size(ins, 1);
        tmp(i) = max(ins(i, :));
    end
    figure;plot(tmp)
    length(find(tmp>1e-1))
end

end

