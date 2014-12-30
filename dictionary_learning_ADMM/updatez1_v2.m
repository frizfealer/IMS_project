function [ uZ1 ] = updatez1_v2( Z1, U1, W, rho, lambda )
%updatez1_v2 update z1 with second ADMM formulation:
%z1 = W
W = W(:);
U1 = U1(:);
C = speye(length(W));
d = W-1/rho*U1;

%Starting point
opts.init=1;      
opts.x0 = Z1(:);
% termination criterion
opts.tFlag=5;       % run  to Xtol
% opts.maxIter=200;   % maximum number of iterations
opts.tol=1e-4;
% normalization
opts.nFlag=0;       % without normalization
opts.rFlag=0;       % the input parameter 'rho' is a ratio in (0, 1)
opts.lFalg = 1;
lambda = lambda/rho;
opts.rsL2 = 0;
[tmp, funVal]=LeastR(C, d, lambda, opts);
uZ1 = reshape( tmp,  size(Z1) );




% options = glmnetSet;
% % options.nlambda=5;
% options.intr=0;
% options.standardize=0;
% options.thresh = 1e-16;
% % options.lambda_min = 0;
% lambda = lambda / rho / size(C, 1);
% options.lambda = [lambda*3, lambda*2, lambda, lambda/2, lambda/3];
% fit = glmnet( C,d, 'gaussian', options );
% tmp = fit.beta(:, 3);
% tmp = reshape( tmp, size(Z1) );
% uZ1 = tmp;
end

