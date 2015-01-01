function [ dfVec ] = testThetaMax( Y, D, startVal, intVal, intThres, testNum )
%testLambdaLimit test lambda max
endVal = startVal+testNum*intVal;
cnt = 1;
for theta = startVal:intVal:endVal
    fprintf( 'theta = %g\n', theta );
    [~, hei, wid] = size(Y);
    aMatrix = ones( hei, wid );
    itNum = 100;
    L1Flag = 1;
    logFY = [];
    initVar = [];
    [WResStruct] = updateW_ADMM_v3( Y, D, aMatrix, itNum, 1e-32, theta, L1Flag, logFY, initVar, 0, 5*1e-2 );
    %figure; subplot(1, 2, 1); plot(WResStruct.LPAry(:, 1)); subplot(1, 2, 2); plot(WResStruct.LPAry(:, 2) );
    BlkDS = conBLKDS( Y );
    ins = WResStruct.W(:, BlkDS.indMap == 1 );
    tmp = zeros( size(ins, 1), 1);
    for i = 1:size(ins, 1);
        tmp(i) = max(ins(i, :));
    end
    figure;plot(tmp)
    dfVec(cnt) = length(find(tmp>intThres));
    cnt = cnt + 1;
    if cnt == 2
        if dfVec(cnt) == dfVec(cnt-1)
            fprintf( 'df is not changed.\n' );
        end
        if dfVec(cnt) == 0
            break;
        end
    end
end

end

