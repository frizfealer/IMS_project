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
    lambda = 1e-32;
    theta = 1e-32;
    L1Flag = 1;
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
    

else
end

end

