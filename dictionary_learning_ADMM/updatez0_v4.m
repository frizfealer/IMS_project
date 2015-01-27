function [ z0 ] = updatez0_v4( alphaFlag, LINK_FUNC, Y, Z0_init, D, W, W0, U0, rho, scaleFac, varargin  )
%updatez0_i,j using minfunc function to implement
%y (location i, j) [s 1], z0 (location i, j) [s 1], D[s m] 
%w (location i, j) [m 1], w0 (location i, j) scalar, u0 (location i, j)[s 1]
%[ val ] = z0_i_j_termTest( y, z0, D, w, w0, u0, rho )

[sLen, nLen] = size(Y);
res1 = - D*W - repmat( W0(:)', sLen, 1 ) + 1/rho*U0;
if alphaFlag == 0
    z0 = -res1;
elseif alphaFlag == 1
    Y = Y(:);
    res1 = res1(:);
    coordXY = 1:sLen*nLen;
    vNum = sLen*nLen;
    if strcmp( LINK_FUNC, 'log' ) == 1
        targetFunc = @(Z0) Z0_termFunc_log( Y, Z0, rho, res1, scaleFac, vNum, coordXY );
    elseif strcmp( LINK_FUNC, 'identity' ) == 1
        targetFunc = @(Z0) Z0_termFunc_identity( Y, Z0, rho, res1, scaleFac, vNum, coordXY );
    elseif strcmp( LINK_FUNC, 'log_gaussain' ) == 1
        targetFunc = @(Z0) Z0_termFunc_log_gaussain( Y, Z0, rho, res1, scaleFac, vNum, coordXY );
    end
    if ~isempty( Z0_init )
        Z0Start = Z0_init(:);
    else
        Z0Start = 1e-6*zeros( sLen*nLen, 1 );
        %     oset = D*w + w0 - 1/rho*u0;
        %     [B,~]= lassoglm( eye(sLen), y,'poisson', 'Offset', oset, 'Alpha', 1e-8, 'Lambda', rho, 'RelTol', 1e-8 );
        %     z0 = B+oset;
        %     z0Start = z0; z0Start(z0Start<0) = 1e-8;
    end
    if length(varargin) == 1
        itNum = varargin{1};
    else
        itNum = 1000;
    end
    if length(varargin) == 2
        options.Method = varargin{2};
    else
        options.Method = 'newton';
    end
%     lowerB = zeros(sLen,1);
%     upperB = zeros(sLen,1); upperB(:) = Inf;
%     options = optimoptions('fminunc');
%     options.MaxFunEvals= 1e6;
%     options.MaxIter = itNum;
%     options.TolFun=1e-6;
%     options.TolX=1e-6;
%     options.Algorithm='trust-region';
    options.Display='none';
%     options.Hessian='on';
% %     options.Diagnostics = 'off';
% %     options.HessUpdate='steepdesc';
%      % options.UseParallel='always';
%     options.GradObj='on';
%     [ z0, ~ ] = fminunc( targetFunc, Z0Start, options );
     [ z0, ~ ] = minFunc( targetFunc, Z0Start, options );
    z0 = reshape( z0, sLen, nLen );
end

end

function [ val,grad, H ] = Z0_termFunc_log( Y, Z0, rho, res1, scaleFac, vNum, coordXY )
%Z0_termFunc the terms consisted of z0 in ADMM formulation
eZ0 = exp(Z0);
%     fprintf( 'll = %g, pp = %g\n', full(scaleFac*sum( ( -y.*z0 + eZ0 ))), full(rho / 2 * sum( (z0 +res1).^2 )));
val = sum( scaleFac*( -Y.*Z0 + eZ0 ) ) + rho / 2 * sum( (Z0 +res1).^2 );
grad = scaleFac*(-Y + eZ0) + rho* ( Z0 + res1 );
tmp = scaleFac*eZ0 + rho;
H = sparse( coordXY, coordXY, tmp, vNum, vNum );
end

function [ val,grad, H ] = Z0_termFunc_identity( Y, Z0, rho, res1, scaleFac, vNum, coordXY )
%Z0_termFunc the terms consisted of z0 in ADMM formulation
%     fprintf( 'll = %g, pp = %g\n', full(scaleFac*sum( ( -y.*z0 + eZ0 ))), full(rho / 2 * sum( (z0 +res1).^2 )));
val = sum( scaleFac*( -Y.*log(Z0) + Z0 ) ) + rho / 2 * sum( (Z0 +res1).^2 );
grad = scaleFac*(-Y./(Z0) + 1) + rho* ( Z0 + res1 );
tmp = scaleFac*(Y./(Z0).^2) + rho;
H = sparse( coordXY, coordXY, tmp, vNum, vNum );
end

function [ val,grad, H ] = Z0_termFunc_log_gaussain( Y, Z0, rho, res1, scaleFac, vNum, coordXY )
residual = log(Y+1e-32) - Z0;
val = sum( scaleFac*(residual).^2 ) + rho/2 * sum( (Z0+res1).^2 );
grad = -2*scaleFac*(residual) + rho*(Z0 + res1);
tmp = 2 + rho;
H = sparse( coordXY, coordXY, tmp, vNum, vNum );
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


