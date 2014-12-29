function [ uZ1 ] = updatez1_v2( Z1, U1, W, rho, lambda )
%updatez1_v2 update z1 with second ADMM formulation:
%z1 = W
W = W(:);
U1 = U1(:);
C = speye(length(W));
d = W-1/rho*U1;
options = glmnetSet;
% options.nlambda=5;
options.intr=0;
options.standardize=0;
options.thresh = 1e-16;
% options.lambda_min = 0;
lambda = lambda / rho / size(C, 1);
options.lambda = [lambda*3, lambda*2, lambda, lambda/2, lambda/3];
fit = glmnet( C,d, 'gaussian', options );
uZ1 = fit.beta(:, 3);
end

