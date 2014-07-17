function [ val, grad] = exponetialPriorFunc( x, Y, B, logAlpha, logLAlpha, thirdT, alpha, beta )
%exponetialPriorFunc exponetialPrior function with linear regression
dataLen = length(x);
sumLogX = sum(log(x));
resY = Y-B*x;
val = 1/2*sum((resY).^2) + dataLen*( -logAlpha - logLAlpha + thirdT ) - (-alpha-1)*sumLogX*beta;
fprintf('%g\n', dataLen*( -logAlpha - logLAlpha + thirdT ) - (-alpha-1)*sumLogX );
grad = -B'*resY - (-alpha-1)./x*beta;
end

