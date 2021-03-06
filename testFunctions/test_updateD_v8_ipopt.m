function [ output_args ] = test_updateD_v8_ipopt( testCaseNum )
%test_updateD_v8_ipopt test function for updateD_v8_ipopt
if testCaseNum == 1
     %% synthesize data firs
    SLEN = 10;
    MLEN = 10;
    HEIGHT = 30;
    WIDTH = 30;
    DOptions.sparsePrec = 1;
    DOptions.coheMax = 1;
    WOptions.supPrec = 1;
    WOptions.sparsePrec = 1;
    WOptions.scale = 5000;
    verbose = 0;
    LINK_FUNC = 'identity';
    CONSTRAINTS = 'L1';
    [simData] = synthesizeData_Poisson( LINK_FUNC, CONSTRAINTS, SLEN, MLEN, HEIGHT, WIDTH, [], DOptions, WOptions, verbose );
    aMatrix = ones( HEIGHT, WIDTH );
    %% itNum = 100
    itNum = 200;
    phi = 1e-32;
    scaleFac = 1;
    [ D ] = initD( simData.gY, simData.uDTemplate, 'NNMF', [] );
    [ preVal, ~, ~ ] = LP_DL_Poiss( LINK_FUNC, aMatrix, simData.gY, simData.gW, simData.gW0, simData.gD, 1e-32, phi, 1e-32, scaleFac, [] );
    [ finalD ] = updateD_v8_ipopt( LINK_FUNC, CONSTRAINTS, simData.gY, simData.gW, simData.gW0, simData.gD, simData.uDTemplate, aMatrix, 0, phi, scaleFac, itNum, 1e-3 );
    [ aftVal, ~, ~ ] = LP_DL_Poiss( LINK_FUNC, aMatrix, simData.gY, simData.gW, simData.gW0, finalD, 1e-32, phi, 1e-32, scaleFac, [] );

    figure;subplot(1,2,1);imagesc(simData.gD);colorbar;subplot(1,2,2);imagesc(finalD);colorbar;
    aftVal-preVal
elseif testCaseNum == 2
     %% synthesize data firs
    SLEN = 10;
    MLEN = 10;
    HEIGHT = 30;
    WIDTH = 30;
    DOptions.sparsePrec = 1;
    DOptions.coheMax = 1;
    WOptions.supPrec = 1;
    WOptions.sparsePrec = 1;
    WOptions.scale = 1;
    verbose = 0;
    LINK_FUNC = 'log_gaussain';
    CONSTRAINTS = 'L1';
    [simData] = synthesizeData_Poisson( LINK_FUNC, CONSTRAINTS, SLEN, MLEN, HEIGHT, WIDTH, [], DOptions, WOptions, verbose );
    aMatrix = ones( HEIGHT, WIDTH );
    %% itNum = 100
    itNum = 200;
    phi = 1e-32;
    scaleFac = 1;
    [ D ] = initD( simData.gY, simData.uDTemplate, 'NNMF', [] );
    [ preVal, ~, ~ ] = LP_DL_Poiss( LINK_FUNC, aMatrix, simData.gY, simData.gW, simData.gW0, simData.gD, 1e-32, phi, 1e-32, scaleFac, [] );
    [ finalD ] = updateD_v8_ipopt( LINK_FUNC, CONSTRAINTS, simData.gY, simData.gW, simData.gW0,D, simData.uDTemplate, aMatrix, 0, phi, scaleFac, itNum, 1e-3, [] );
    [ aftVal, ~, ~ ] = LP_DL_Poiss( LINK_FUNC, aMatrix, simData.gY, simData.gW, simData.gW0, finalD, 1e-32, phi, 1e-32, scaleFac, [] );

    figure;subplot(1,2,1);imagesc(simData.gD);colorbar;subplot(1,2,2);imagesc(finalD);colorbar;
    aftVal-preVal

elseif testCaseNum == 3
     %% synthesize data firs
    SLEN = 10;
    MLEN = 10;
    HEIGHT = 30;
    WIDTH = 30;
    DOptions.sparsePrec = 1;
    DOptions.coheMax = 1;
    WOptions.supPrec = 1;
    WOptions.sparsePrec = 1;
    WOptions.scale = 1;
    verbose = 0;
    LINK_FUNC = 'negative_binomial';
    CONSTRAINTS = 'L1';
    [simData] = synthesizeData_Poisson( LINK_FUNC, CONSTRAINTS, SLEN, MLEN, HEIGHT, WIDTH, [], DOptions, WOptions, verbose );
    aMatrix = ones( HEIGHT, WIDTH );
    %% itNum = 100
    itNum = 200;
    phi = 1e-32;
    scaleFac = 1;
    [ D ] = initD( simData.gY, simData.uDTemplate, 'NNMF', [] );
    [ preVal, ~, ~ ] = LP_DL_Poiss( LINK_FUNC, aMatrix, simData.gY, simData.gW, simData.gW0, simData.gD, 1e-32, phi, 1e-32, scaleFac, [], [], 1e-2 );
    [ finalD, kappa ] = updateD_v8_ipopt( LINK_FUNC, CONSTRAINTS, simData.gY, simData.gW, simData.gW0, D, simData.uDTemplate, aMatrix, 0, phi, scaleFac, itNum, 1e-3, [] );
    [ aftVal, ~, ~ ] = LP_DL_Poiss( LINK_FUNC, aMatrix, simData.gY, simData.gW, simData.gW0, finalD, 1e-32, phi, 1e-32, scaleFac, [], [], kappa );

    figure;subplot(1,2,1);imagesc(simData.gD);colorbar;subplot(1,2,2);imagesc(finalD);colorbar;
    aftVal-preVal

end

end

