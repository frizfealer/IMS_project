function [ D_NNMF, Res ] = NNMF_als_wrapper( inY, initD, BlkDS, mNum, rNum, paraFlag, implementedMethod )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
opt = statset( 'MaxIter', 200, 'Display', 'final', 'TolFun', 1e-16, 'TolX', 1e-16, 'UseParallel', paraFlag );
algoUsed = 'als';
if ~isempty(initD)
    initDFlag = true;
    assert( mNum == size(initD, 2) );
else
    initDFlag = false;
end


if strcmp( implementedMethod, 'default' ) == 1
    if initDFlag == true
        [W,H, Res] = nnmf( inY(:, BlkDS.indMap==1), mNum, 'w0', initD, 'replicates', rNum,...
            'options', opt,...
            'algorithm', algoUsed );
    else
        [W,H, Res] = nnmf( inY(:, BlkDS.indMap==1), mNum, 'replicates', rNum,...
            'options', opt,...
            'algorithm', algoUsed );
    end
elseif strcmp( implementedMethod, 'NMF_MATLAB_TOOLBOX' ) == 1
     option.residual = 1e-16;
     option.tof = 1e-16;
     [ W, H, numIter, tElapsed, finalResidual]=nmfnnls( inY(:, BlkDS.indMap==1 ), mNum, option );
end
for i = 1:mNum
    W(:, i) = W(:, i) / norm( W(:, i) );
end
D_NNMF = W;

end

