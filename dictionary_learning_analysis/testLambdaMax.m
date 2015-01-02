function [ dfVec, lambdaVec ] = testLambdaMax( Y, D, startVal, intVal, intThres, testNum, seq, w_tol )
%testLambdaMax test lambda max
if ~isempty(seq)
    lambdaVec = seq;
else
    endVal = startVal+testNum*intVal;
    lambdaVec = startVal:intVal:endVal;
end
dfVec = zeros(length(lambdaVec), 1);
for i = 1:length(lambdaVec)
    lambda = lambdaVec(i);
    fprintf( 'lambda = %g\n', lambda );
    [~, hei, wid] = size(Y);
    aMatrix = ones( hei, wid );
    itNum = 100;
    theta = 1e-32;
    L1Flag = 1;
    logFY = [];
    initVar = [];
    scaleFactor = computeScaleFactor( Y, aMatrix );
    [WResStruct] = updateW_ADMM_v3( Y, D, aMatrix, itNum, lambda, theta, L1Flag, logFY, initVar, scaleFactor, 0, w_tol );
    BlkDS = conBLKDS( Y );
    %figure; subplot(1, 2, 1); plot(WResStruct.LPAry(:, 1)); subplot(1, 2, 2); plot(WResStruct.LPAry(:, 2) );
    ins = WResStruct.W(:, BlkDS.indMap == 1 );
    tmp = zeros( size(ins, 1), 1);
    for j = 1:size(ins, 1);
        tmp(j) = max(ins(j, :));
    end
    %figure;plot(tmp)
    dfVec(i) = length(find(tmp>intThres));
    %save( 'testLambdaMax_temp.mat', 'dfVec', 'lambdaVec' );
end

end

