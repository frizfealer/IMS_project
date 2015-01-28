function [ maxLambda, minLambda ] = estimateLambdaMaxMin( Y, D, W, W0, linkFunc )
%estimateLambdaMax 
[sLen, ~, ~]= size(Y);
preY = D*W(:,:)+repmat( W0(:)', sLen, 1 );
if strcmp( linkFunc, 'identity' ) == 1
    val1 = sum(Y(:).*log(preY(:)+1e-32)-preY(:) );
    val2 = sum(W(:));
    maxLambda = val1/ val2*10;
    tmp = 1e-4*( val1/ (val2*1*10^-(ceil(log10(val2))) ) );
    minLambda = 1*10^-(ceil(log10(val2)))*tmp;
elseif strcmp( linkFunc, 'log' ) == 1
    val1 = sum( Y(:).*preY(:)-exp(preY(:)) );
    val2 = sum(W(:));
    maxLambda = val1/val2*10;
    tmp = 1e-4*( val1/ (val2*1*10^-(ceil(log10(val2))) ) );
    minLambda = 1*10^-(ceil(log10(val2)))*tmp;
elseif strcmp( linkFunc, 'log_gaussain' ) == 1
    val1 = sum((log(Y(:)+1e-32)-preY(:)).^2);
    val2 = sum(W(:));
    maxLambda = val1/val2*10;
    tmp = 1e-4*( val1/ (val2*1*10^-(ceil(log10(val2))) ) );
    minLambda = 1*10^-(ceil(log10(val2)))*tmp;
end


end

