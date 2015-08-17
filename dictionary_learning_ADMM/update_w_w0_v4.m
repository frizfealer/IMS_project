function [ w_w0 ] = update_w_w0_v4( z0, u0, z1, u1, z2, u2, CforUW, rho, w_init, offsetFlag )
%update_w_w0_ver3 update w and w0 in the whole sample
%11172014 write up
%using glmnet, the speed is at 0.7s, with five lambda value
%using lsqlin, the speed is at 16s
%using SLEP, the speed is at 12s

w_init = w_init(:);
% mLen = length(w_init);
nLen = size(z0, 2);
if offsetFlag == 1
    z2 = [z2 zeros( size(z2, 1), 1 ) ];
    u2 = [u2 zeros( size(u2, 1), 1 ) ];
    tmpD = [z0+1/rho*u0; z1+1/rho*u1; sparse(1, nLen) ];
else
    z2 = [z2 ];
    u2 = [u2];
    tmpD = [z0+1/rho*u0; z1+1/rho*u1; ];
end

d = tmpD(:);
tmpD = (z2+1/rho*u2)';
d = [d; tmpD(:)];
% options = optimoptions('lsqlin');
% options.MaxIter = 200;
% options.TolFun=1e-8;
% options.Display='final';
% lb = zeros(mLen, 1);
% tic
% [w_w0,~,~,~, ~] = lsqlin( C, d, [], [], [], [], lb, [], w_init, options );
% toc
% 
% opts=[];
% 
% % % Starting point
opts.init=1;      
opts.x0 = w_init;
% termination criterion
opts.tFlag=5;       % run  to Xtol
opts.maxIter=200;   % maximum number of iterations
opts.tol=1e-4;
% normalization
opts.nFlag=0;       % without normalization
opts.rFlag=1;       % the input parameter 'rho' is a ratio in (0, 1)
opts.rsL2=0;        % the squared two norm term
opts.lFalg = 1;
% C = sparse(C);
% tic
% [w_w0, funVal]=nnLeastR(CforUW, d, 0, opts);
% toc
% fprintf( '%d\n', norm(x-w_w0, 2 ) );
% % w_w0 = x;
options = glmnetSet;
options.cl=[0;inf];
% options.nlambda=5;
options.lambda = [3 2 1 0];
options.intr=0;
options.standardize=0;
%options.thresh = 1e-16;
options.thresh = 1e-8;
options.lambda_min = 0;
options.ltype='modified.Newton'; %is faster then 'Newton'
fit = glmnet( CforUW,d, 'gaussian', options );
w_w0 = fit.beta(:, end);

end

function C = replicateC( inD, repNum )
[i, j] = find( inD ~= 0 );
[sLen, mLen] = size(inD);
valVec = nonzeros( inD );
valVec = repmat( valVec, repNum, 1 );
yVec = zeros(length(i)*repNum, 1);
xVec = zeros(length(i)*repNum, 1);
unitLen = length(i);
for it = 1:repNum
    xVec( ((it-1)*unitLen+1):(unitLen*it) ) = j + ((it-1)*mLen);
    yVec( ((it-1)*unitLen+1):(unitLen*it) ) = i + ((it-1)*sLen);
end
C = sparse( yVec, xVec, valVec, sLen*repNum, mLen*repNum );

% for it = 2:repNum
%     tmp = S( ((it-1)*sLen+1):(it*sLen), : ) - S(1:sLen, :);
%     if ~isempty( find( tmp ~= 0, 1 ) )
%         fprintf( '%d\n', it );
%     end
% end
end
