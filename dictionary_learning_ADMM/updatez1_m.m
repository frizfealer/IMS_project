function [ z1m, warmZ1m ] = updatez1_m( hei, wid, rho, Wm, u1m, lambda, theta, D, z1Start )
%updatez1_m update z1m (that is the m element across all grids of z1)
%rho: scalar
%Wm: outW(m,:) size( 1, hei*wid), u1m: u1(m,:), same size as Wm
%according to 02062014 writeup

% %build up E matrix, the relationship between the Wm in the grid
% %The |w_i,j-w_i-1,j| relationship
% interval = hei*wid;
% yAry = zeros( 2, wid*(hei-1) );
% xAry = zeros( 2, wid*(hei-1) );
% dataAry = ones( 2, wid*(hei-1) );
% dataAry(1,:) = -1;
% cnt = 1;
% for i=1:(interval-1)
%     if mod( i, hei ) ~= 0
%         yAry(:, cnt) = i;
%         xAry(1, cnt) = i;
%         xAry(2, cnt) = i + 1;
%         cnt = cnt + 1;
%     end
% end
% yAry = yAry(:);
% xAry = xAry(:);
% dataAry = dataAry(:);
% E1 = sparse( yAry, xAry, dataAry, interval, interval );
% %The |w_i,j-w_i,j-1| relationship
% interval = hei*wid;
% yAry = zeros( 2, hei*(wid-1) );
% xAry = zeros( 2, hei*(wid-1) );
% dataAry = ones( 2, hei*(wid-1) );
% dataAry(1,:) = -1;
% for i = 1:(interval-hei)
%     yAry(:,i) = i;
%     xAry(1,i) = i;
%     xAry(2,i) = i+hei;
% end
% yAry = yAry(:);
% xAry = xAry(:);
% dataAry = dataAry(:);
% E2 = sparse( yAry, xAry, dataAry, interval, interval );

% tic;
rrho = rho ^0.5;
y = (rrho*Wm - 1/rrho*u1m)'; % a hei*wid-by-1 vector
X = sparse( 1:length(y), 1:length(y), rrho );
% D = [E1; E2];
lAry = ones( size(X, 1), 1 )*lambda+1e-8;
tAry = ones( size(D, 1), 1 )*theta+1e-8;
% toc;
[w,warm,gap, it1, it2, condM]=lipm2m_ADMM( y, X, D, lAry, tAry, z1Start, 0 );
% fprintf( 'condM: %g, gap: %g it1: %d it2: %d\n', condM, gap, it1, it2 );
z1m = w';
warmZ1m = warm;

% options = optimoptions('fminunc');
% options.MaxFunEvals= 1e6;
% options.MaxIter = 200;
% options.TolFun=1e-6;
% options.TolX=1e-6;
% options.Algorithm='quasi-newton';
% options.Display='none';
% options.Diagnostics = 'off';
% test = @(W) z1_m_termFunc( y, W, D, lambda, theta );
% [w,~] = fminunc(test,z1Start, options);
% z1m = w';

end

function [ val ] = z1_m_termFunc( y,W, D, l1, l2 )
%(1/2)||y - Xw||^2_2 + l1'||w||_1 + l2'||Dw||_1
val = norm(y(:)-W(:)) + l1*abs(sum(W(:))) + l2*abs(sum(sum(D*W')));
end