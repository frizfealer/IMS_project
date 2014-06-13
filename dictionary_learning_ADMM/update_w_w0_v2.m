function [ w_w0 ] = update_w_w0_v2( z0, u0, z1, u1, CforUW, rho, w_init )
%update_w_w0_ver2 update w and w0 at each location in
%02062014 write up

mLen = size( z1, 1 );
d = [z0+1/rho*u0; z1+1/rho*u1; 0];
C = CforUW;
options = optimoptions('lsqlin');
options.MaxIter = 100;
options.TolFun=1e-6;
options.Display='none';
lb = zeros(mLen+1, 1);
[w_w0,~,~,~, ~] = lsqlin( C, d, [], [], [], [], lb, [], w_init, options );
% 
opts=[];

% % Starting point
% opts.init=1;      
% opts.x0 = w_init;
% % termination criterion
% opts.tFlag=10;       % run .maxIter iterations
% opts.maxIter=100;   % maximum number of iterations
% opt.tol=1e-5;
% % normalization
% opts.nFlag=0;       % without normalization
% opts.rFlag=1;       % the input parameter 'rho' is a ratio in (0, 1)
% opts.rsL2=0;        % the squared two norm term
% 
% C = sparse(C);
% tic
% [x, funVal]=nnLeastR(C, d, 1e-6, opts);
% toc
% fprintf( '%d\n', norm(x-w_w0, 2 ) );
% w_w0 = x;


end

