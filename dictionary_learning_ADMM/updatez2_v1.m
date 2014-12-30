function [ uZ2 ] = updatez2_v1( Z2, U2, rho, R, W, theta, L1Flag )
%updatez2_v1 update z2 in the second ADMM formulation:
%z2 = R*W';
diffW = R*W';
diffW = diffW'; diffW = diffW(:);
U2 = U2'; U2 = U2(:);
C = speye(length(diffW));
d = diffW-1/rho*U2;
%Starting point
opts.init=1;      
opts.x0 = Z2(:);
% termination criterion
opts.tFlag=5;       % run  to Xtol
% opts.maxIter=200;   % maximum number of iterations
opts.tol=1e-4;
% normalization
opts.nFlag=0;       % without normalization
opts.rFlag=0;       % the input parameter 'rho' is a ratio in (0, 1)
opts.lFalg = 1;
if L1Flag == 1
    lambda = theta/rho;
    opts.rsL2 = 0;
elseif L1Flag == 0
    lambda = 0;
    opts.rsL2 = theta*2/rho;      % the squared two norm term
end
[tmp, funVal]=LeastR(C, d, lambda, opts);
tmp = reshape( tmp,  size(W, 1), size(R, 1) );
uZ2 = tmp';

% options = glmnetSet;
% % options.nlambda=5;
% options.intr=0;
% options.standardize=0;
% options.thresh = 1e-32;
% % options.lambda_min = 0;
% oTheta = theta;
% theta = theta / rho / size(C, 1);
% if L1Flag == 0
% %     options.alpha = 0;
% %     theta = theta * 2;
% %     options.lambda_min = theta;
%     theta = 1e-32;
%     E = speye(length(diffW))*sqrt(oTheta);
%     C = rowCombSparseMatrix( C, E );
%     d = [d; zeros(length(diffW), 1)];
% end
% options.lambda = [theta*5, theta*4, theta*3, theta*2, theta];
% tic; fit = glmnet( C,d, 'gaussian', options ); toc;
% tmp = fit.beta(:, end);
% tmp = reshape( tmp,  size(W, 1), size(R, 1) );
% uZ2 = tmp';

end

function C = rowCombSparseMatrix( A, B )
    [i, j] = find( A ~= 0 );
    [sLenA, mLenA] = size(A);
    [q, r] = find( B ~= 0 );
    [sLenB, mLenB] = size(B);
    if mLenA ~= mLenB
        fprintf( 'row combine sparse matrix error.\n' );
        C = [];
    else
        valVecA = nonzeros( A );
        valVecB = nonzeros( B );
        allVec = [valVecA; valVecB];
        q = q + sLenA;
        allYAxis = [i; q];
        allXAxis = [j; r];
        C = sparse( allYAxis, allXAxis, allVec, sLenA+sLenB, mLenA );
    end
end

