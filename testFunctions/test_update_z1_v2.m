function [ output_args ] = test_update_z1_v2( testCaseNum )
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
    lambda = 1;
    preZ1 = simData.gW(:, :);
    U1 = zeros( size( preZ1 ) );
    W = simData.gW(:, :);
    fVal(1) = Z1_term_func( W, preZ1, U1, rho, lambda );
    [ uZ1 ] = updatez1_v2( preZ1, U1, W, rho, lambda );
    tmp = reshape( uZ1, MLEN, HEIGHT*WIDTH );
    Z1 = tmp;
    fVal(2) = Z1_term_func( W, Z1, U1, rho, lambda );
    fVal
    fVal = repmat(fVal', 1, 100);
    dfVal = zeros( 2, 100 );
    for lambda = 1:100
        [ Z1 ] = updatez1_v2( preZ1, U1, W, rho, lambda );
        fVal(1, lambda) = Z1_term_func( W, preZ1, U1, rho, lambda );
        fVal(2, lambda) = Z1_term_func( W, Z1, U1, rho, lambda );
        dfVal(1, lambda) = length( find( abs( preZ1(:) ) >= 1e-3 ) );
        dfVal(2, lambda) = length( find( abs( Z1(:) ) >= 1e-3 ) );
    end
    figure; plot(fVal');
    figure; plot(dfVal');
elseif testCaseNum == 2
    load( 'D:\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LN\100831_348_136_12_40hr_0_1XLB_LN_input.mat' );
    load( 'D:\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LN\exp_100831_348_136_12_40hr_0_1XLB_LN_4_res.mat' );
    inY = dataCube(:, :);
    [sLen, hei, wid] = size(dataCube);
    Z1 = expRec.z1;
    U1 = zeros( size( Z1 ) );
    curRho = 1;
    W = expRec.outW;
    lambda = 1e-2;
    fVal = zeros(4, 100);
    dfVal = zeros( 2, 100 );
    for i = 1:100
        for j = 1:BlkDS.blkNum
            loc = BlkDS.B2GMap{j};
            curW = W(:, loc);
            curZ1 = Z1(:, loc);
            preZ1 = curZ1;
            curU1 = U1(:, loc);
            [ curZ1 ] = updatez1_v2( preZ1, curU1, curW, curRho, lambda );
            fVal((j-1)*2+1, i) = Z1_term_func( curW, preZ1, curU1, curRho, lambda );
            fVal((j-1)*2+2, i) = Z1_term_func( curW, curZ1, curU1, curRho, lambda );
            dfVal((j-1)*2+1, i) = length( find( abs(preZ1(:)) >= 1e-3 ) );
            dfVal((j-1)*2+2, i) = length( find( abs(curZ1(:)) >= 1e-3 ) );
        end
        lambda = lambda*(i+1);
    end
    figure; plot(fVal(1:2,:)');
    figure; plot(fVal(3:4,:)');
end

end

function [ val ] = Z1_term_func( W, Z1, U1, rho, lambda )
%this function should be minimized.
val1 = W-1/rho*U1-Z1;
val1 = (norm(val1(:))^2)*rho/2;
val2 = norm( Z1(:), 1 )*lambda;
val = val1+val2;
end

