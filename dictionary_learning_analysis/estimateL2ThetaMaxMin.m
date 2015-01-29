function [ maxTheta, minTheta ] = estimateL2ThetaMaxMin( Y, D, W, W0, linkFunc, varargin )
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
    ttt = log(Y);
    ttt(ttt==-inf) = 0;
    val1 = sum((ttt(:)-preY(:)).^2);
    val2 = (Rall*W(:,:)').^2;
    val2 = sum(val2(:));
    maxTheta = val1/val2*10;
    tmp = 1e-4*( val1/ (val2*1*10^-(ceil(log10(val2))) ) );
    minTheta = 1*10^-(ceil(log10(val2)))*tmp;
elseif strcmp( linkFunc, 'negative_binomial' ) == 1
    kappa = varargin{1};
    val1 = -Y(:)*log(kappa) - Y(:).*preY(:) ...
        + (Y(:)+1/kappa).*log( 1+kappa*exp(preY(:)) ) ...
        - gammaln( Y(:)+1/kappa ) + gammaln( 1/kappa );
    val1 = abs(sum(val1(:)) );
    val2 = (Rall*W(:,:)').^2;
    val2 = sum(val2(:));
    maxTheta = val1/val2*10;
    tmp = 1e-4*( val1/ (val2*1*10^-(ceil(log10(val2))) ) );
    minTheta = 1*10^-(ceil(log10(val2)))*tmp;
end


end

