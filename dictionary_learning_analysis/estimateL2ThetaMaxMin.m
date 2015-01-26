function [ maxTheta, minTheta ] = estimateL2ThetaMaxMin( Y, D, W, W0, linkFunc )
%estimateLambdaMax 
[sLen, hei, wid]= size(Y);
BlkDS = conBLKDS( Y );
% nLen = hei*wid;
[Rall] = genSparseGroupingMatrix( hei, wid, 1 );
preY = D*W(:,:)+repmat( W0(:)', sLen, 1 );
if strcmp( linkFunc, 'identity' ) == 1
    val1 = sum( Y(:).*log((preY(:))+1e-32)-preY(:) );
    val2 = (Rall*W(:,:)').^2;
    val2 = sum(val2(:));
    maxTheta = val1/val2*10;
    tmp = 1e-4*( val1/ (val2*1*10^-(ceil(log10(val2))) ) );
    minTheta = 1*10^-(ceil(log10(val2)))*tmp; 
elseif strcmp( linkFunc, 'log' ) == 1
    val1 = sum( Y(:).*preY(:)-exp(preY(:)) );
    val2 = (Rall*W(:,:)').^2;
    val2 = sum(val2(:));
    maxTheta = val1/val2*10;
    tmp = 1e-4*( val1/ (val2*1*10^-(ceil(log10(val2))) ) );
    minTheta = 1*10^-(ceil(log10(val2)))*tmp;
elseif strcmp( linkFunc, 'log_gaussain' ) == 1
end


end

