function [ z0 ] = updatez0_i_j_3( alpha, y, z0_init, D, w, w0, u0, rho, scaleFac, varargin  )
%updatez0_i,j using lassoglm and fmincon function to implement
%y (location i, j) [s 1], z0 (location i, j) [s 1], D[s m] 
%w (location i, j) [m 1], w0 (location i, j) scalar, u0 (location i, j)[s 1]
%[ val ] = z0_i_j_termTest( y, z0, D, w, w0, u0, rho )


res1 = - D*w - w0 + 1/rho*u0;
if alpha == 0
    z0 = -res1;
elseif alpha == 1
    test = @(z0) z0_i_j_termFunc( y, z0, rho, res1, scaleFac );
    [sLen, ~] = size(D);
    if ~isempty( z0_init )
        z0Start = z0_init;
    else
        z0Start = 1e-6*zeros(sLen,1);
        %     oset = D*w + w0 - 1/rho*u0;
        %     [B,~]= lassoglm( eye(sLen), y,'poisson', 'Offset', oset, 'Alpha', 1e-8, 'Lambda', rho, 'RelTol', 1e-8 );
        %     z0 = B+oset;
        %     z0Start = z0; z0Start(z0Start<0) = 1e-8;
    end
    if length(varargin) == 1
        itNum = varargin{1};
    else
        itNum = 400;
    end
    
%     lowerB = zeros(sLen,1);
%     upperB = zeros(sLen,1); upperB(:) = Inf;
    options = optimoptions('fminunc');
    options.MaxFunEvals= 1e6;
    options.MaxIter = itNum;
    options.TolFun=1e-6;
    options.TolX=1e-6;
    options.Algorithm='trust-region';
    options.Display='none';
    options.Hessian='on';
    options.Diagnostics = 'off';
%     options.HessUpdate='steepdesc';
     % options.UseParallel='always';
    options.GradObj='on';
    [z0,~] = fminunc(test,z0Start, options);
end

end

function [ val,grad, H ] = z0_i_j_termFunc( y, z0, rho,res1, scaleFac )
%z0_i_j_termTest the AL with z0 terms
%y (location i, j) [s 1], z0 (location i, j) [s 1], D[s m] 
%w (location i, j) [m 1], w0 (location i, j) scalar, u0 (location i, j)[s 1]
eZ0 = exp(z0);
%     fprintf( 'll = %g, pp = %g\n', full(scaleFac*sum( ( -y.*z0 + eZ0 ))), full(rho / 2 * sum( (z0 +res1).^2 )));
val = sum( scaleFac*( -y.*z0 + eZ0 ) ) + rho / 2 * sum( (z0 +res1).^2 );
grad = scaleFac*(-y + eZ0) + rho* ( z0 + res1 );
sLen = size(z0, 1);
tmp = scaleFac*eZ0 + rho;
idx = 1:sLen;
H = sparse( idx, idx, tmp, sLen, sLen );
end

function [ hMat ] = z0_i_j_HessianFFunc( z0, lambda, rho  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
sLen = size(z0, 1);
eZ0 = exp(z0);
tmp = eZ0 + rho;
hTmp = zeros( sLen, 1);
for i = 1:sLen
    hTmp(i) = tmp(i);
end
idx = 1:sLen;
hMat = sparse( idx, idx, hTmp, sLen, sLen );
% for i = 1:sLen
%     hMat(i,i) = tmp(i);
% end
% hMat = sparse(hMat);
end
