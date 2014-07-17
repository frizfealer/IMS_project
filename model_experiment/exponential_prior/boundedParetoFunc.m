function [ val, grad ] = boundedParetoFunc( alpha, X, L, H, logL, sumLogX )
%boundedParetoFunc boundedPareto distribution defined in Wiki
%   alpha*L^alpha*x^(-alpha-1)/(1-(L/H)^alpha)
dataLen = length(X);
ratio = L/H;
val = dataLen*( -log(alpha)-logL*alpha + log(1-ratio^alpha) ) - (-alpha-1)*sumLogX;
grad = dataLen*( -1/alpha - logL + 1/(1-L^alpha)*(-logL*L^alpha) ) + sumLogX;
end

