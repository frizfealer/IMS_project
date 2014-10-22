function [ val ] = computeBIC( mode, Y, D, W, W0, logFY )
%computeBIC compute Bayesian information criterion
%logFY: a vector for log factorial Y
%if mode == 1, using the formula -2*lnL+k*ln(n)
%if mode == 2, using the formula n*ln(Y-D*W) + k*ln(n)
%where L is the likelihood value
%k is # unique nonzero elements of features
%n is # data points in x, in our case, is Spectral Size*Sample Size
Y = Y(:, :);
threshold1 = 1e-3;
dfW = length(find(W>threshold1));
threshold2 = 1e-3;
dfD = length(find(D>threshold2));
n = length( Y(:) );
sLen = size( Y, 1 );
preY = D* W(:, :) + repmat( W0(:)', sLen, 1 );
if mode == 1
    val = -2*( sum(sum( Y.*preY - exp( preY ) )) - sum(logFY) ) + dfW*dfD*log(n);
elseif mode == 2
    y = log(Y);
    y(y ==-inf) = 0;
    val = n*log( norm( y-preY, 'fro' )^2/n ) + dfW*dfD*log(n);
end

end

