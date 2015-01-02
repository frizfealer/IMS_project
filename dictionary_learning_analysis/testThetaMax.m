function [ dfVec, thetaVec ] = testThetaMax( Y, D, startVal, intVal, intThres, testNum, seq, w_tol )
%testThetaMax test lambda max
if ~isempty(seq)
    thetaVec = seq;
else
    endVal = startVal+testNum*intVal;
    thetaVec = startVal:intVal:endVal;
end
dfVec = zeros(length(thetaVec), 1);
[~, hei, wid] = size(Y);
BlkDS = conBLKDS( Y );
[Rall] = genSparseGroupingMatrix( hei, wid, 1 );
for i = 1:BlkDS.blkNum
    Rblk{i} = [];
    Rb = Rall(:, BlkDS.B2GMap{i});
    ins = sum(abs(Rb), 2);
    Rb(ins<2,:) = [];
    Rblk{i} = Rb;
end
for i = 1:length(thetaVec)
    theta = thetaVec(i);
    fprintf( 'theta = %g\n', theta );
    aMatrix = ones( hei, wid );
    itNum = 100;
    L1Flag = 1;
    logFY = [];
    initVar = [];
    scaleFactor = computeScaleFactor( Y, aMatrix );
    [WResStruct] = updateW_ADMM_v3( Y, D, aMatrix, itNum, 1e-32, theta, L1Flag, logFY, initVar, scaleFactor, 0, w_tol );
    %figure; subplot(1, 2, 1); plot(WResStruct.LPAry(:, 1)); subplot(1, 2, 2); plot(WResStruct.LPAry(:, 2) );
    ins = Rall*WResStruct.W'; ins = abs(ins);
    tmp = zeros( size(ins, 1), 1);
    for j = 1:size(ins, 1);
        tmp(j) = max(ins(j, :));
    end
    %figure;plot(tmp)
    dfVec(i) = length(find(tmp>intThres));
end

end

