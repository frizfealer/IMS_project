function [ devVec, lambdaVec ] = testLambdaMaxL2( Y, D, startVal, intVal, testNum, seq, w_tol )
%testLambdaMax test lambda max
if ~isempty(seq)
    lambdaVec = seq;
else
    endVal = startVal+testNum*intVal;
    lambdaVec = startVal:intVal:endVal;
end
devVec = zeros(length(lambdaVec), 1);
for i = 1:length(lambdaVec)
    lambda = lambdaVec(i);
    fprintf( 'lambda = %g\n', lambda );
    [~, hei, wid] = size(Y);
    aMatrix = ones( hei, wid );
    itNum = 100;
    theta = 1e-32;
    L1Flag = 0;
    logFY = [];
    initVar = [];
    scaleFactor = computeScaleFactor( Y, aMatrix );
    [WResStruct] = updateW_ADMM_v3( Y, D, aMatrix, itNum, lambda, theta, L1Flag, logFY, initVar, scaleFactor, 0, w_tol );
    BlkDS = conBLKDS( Y );
    %figure; subplot(1, 2, 1); plot(WResStruct.LPAry(:, 1)); subplot(1, 2, 2); plot(WResStruct.LPAry(:, 2) );
    ins = WResStruct.W(:, BlkDS.indMap == 1 );
    tmp = zeros( size(ins, 1), 1);
    for j = 1:size(ins, 1);
        tmp(j) = mean(ins(j, :));
    end
    %figure;plot(tmp)
    thres = 0.66*length(tmp);
    [counts, centers] = hist(tmp, 1000);
    accCnt = 0;
    for j = 1:1000
        accCnt = accCnt + counts(j);
        if accCnt >= thres
            break;
        end
    end
    devVec(i) = centers(j);
    %figure; hist(tmp, 1000);
%     dfVec(i) = length(find(tmp>intThres));
    %save( 'testLambdaMax_temp.mat', 'dfVec', 'lambdaVec' );
end

end

