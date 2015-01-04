function [WResStruct] = updateW_ADMM_testing( Y, D, itNum, logFY, initVar, scaleFactor, varargin )

%% variable setting

if ~isempty(varargin)  
    if length(varargin) == 2
        Wtol = varargin{2};
    else
        Wtol = 1e-3;
    end
else
    Wtol = 1e-3;
end

fprintf( 'iteration number for W_update_ADMM = %d\n', itNum );
[sLen, mLen] = size(D);
[~, nLen] = size( Y(:, :) );
%% initializing z0 as log(inY)

if isempty( initVar )
    z0 = log(Y); z0(z0==-inf)=0; %initialize zo as log(Y)
    W = zeros( mLen, nLen ); W0 = zeros( 1, nLen );
else
    z0 = initVar.z0(:, :);
    W = initVar.W(:, :); W0 = initVar.W0(:, :);
end

LPAryADMM = zeros( itNum+1, 1 );
resRecAry = zeros( itNum, 4, 1 );
curRhoAry = zeros( itNum, 1 );
curRhoAry(1, :) = 1;
u0 = sparse( zeros( size( z0(:, :) ) ) );

moutD = [D ones(sLen, 1)];
aPartofCforUW = sparse( [moutD] );
CforUW = replicateC( aPartofCforUW, nLen );
%% main optimization process
    LPAryADMM(1) = LLfunc( Y, W, W0, D, scaleFactor);
    for itNumADMM = 1:itNum
%         tic;
        preW = W;
        fprintf( 'LP: %g ', LPAryADMM(itNumADMM) );
        curRho = curRhoAry(itNumADMM);
        %% update W
        fprintf( 'updating W... ' );
        tmpD = z0+1/curRho*u0;
        d = tmpD(:);
        options = glmnetSet;
        options.cl=[0;inf];
        % options.nlambda=5;
        options.lambda = [3 2 1 0];
        options.intr=0;
        options.standardize=0;
        %options.thresh = 1e-16;
        options.thresh = 1e-8;
        options.lambda_min = 0;
        fit = glmnet( CforUW,d, 'gaussian', options );
        w_w0 = fit.beta(:, end);
        tmp = reshape(w_w0, mLen+1, nLen);
        W = tmp(1:(end-1), :);
        W0 = tmp(end,:);
        %% update Z0
        fprintf( 'updating z0... ' );
        preZ0 = z0;
        tmpY = Y(:);
        res1 = - D*W - repmat( W0, sLen, 1 ) + 1/curRho*u0;
        res1 = res1(:);
        coordXY = 1:sLen*nLen;
        vNum = sLen*nLen;
        targetFunc = @(Z0) Z0_termFunc( tmpY, Z0, curRho, res1, scaleFactor, vNum, coordXY );
        Z0Start = z0(:);
        options.Method = 'newton';
        options.Display = 'off';
        [ z0, ~ ] = minFunc( targetFunc, Z0Start, options );
        z0 = reshape( z0, sLen, nLen );
        %% update u0, u1
        fprintf( 'updating u0... ' );
        preY = D*W(:,:) + repmat( W0(:)', sLen, 1);
        %primal residual r0
        r0 = z0 - preY;
        u0 = u0 + curRho*r0;
        %compute the primal residual all
        rf = norm( r0(:) );
        pDim = length( r0(:) );
        %compute the dual variable to see if it is converge
        %dual residual s0
        s0 = curRho*( D'*( z0 - preZ0 ) );
        sf = norm( s0(:) );
        nDim = length( s0(:) );
        
        EPS_ABS = 1e-4;
        EPS_REL = 1e-3;
        
        relEpsPri = max( [norm( z0(:) ) norm( preY(:) )] );
        epsPri = sqrt(pDim)*EPS_ABS + relEpsPri*EPS_REL;
        
        tmp = D'*u0;
        tmp = [tmp; sum(u0)];
        maxRel = norm( tmp(:) );
        epsDual = sqrt(nDim)*EPS_ABS + maxRel*EPS_REL;
        
        resRecAry(itNumADMM, :) = [rf, epsPri, sf, epsDual];
        tmp = full( max( abs( W(:) -  preW(:) ) ) );
        fprintf( '%d: %g %g %g %g %g %g \n',itNumADMM, rf, epsPri, sf, epsDual, curRhoAry(itNumADMM), full(tmp) );
%         time = toc;
%         fprintf( 'time: %g\n', time );        
        %% breaking condition
        if ( rf < epsPri && sf < epsDual ) || ( tmp < Wtol )
                %( max( abs(rf(:)) ) < 1e-3 && max( abs(sf(:)) ) < 1e-3 ) || ...
            break;
        end
        %% update the next stage rho
        if rf > 10*sf
            curRhoAry(itNumADMM+1) = 2 * curRhoAry(itNumADMM);
        elseif sf > 10*rf
            curRhoAry(itNumADMM+1) = 0.5 * curRhoAry(itNumADMM);
        else
            curRhoAry(itNumADMM+1) = curRhoAry(itNumADMM);
        end
        
        LPAryADMM(itNumADMM+1) = LLfunc( Y, W, W0, D, scaleFactor);
    end
    WResStruct.W = sparse( W );
    WResStruct.W0 = sparse( W0 );
    WResStruct.z0 = sparse( z0(:, :) );
    WResStruct.u0 = sparse(u0(:, :) );
    WResStruct.LPAry = LPAryADMM;
    WResStruct.rhoAry = curRhoAry;
    WResStruct.resAry = resRecAry;
    WResStruct.WDiff = tmp;
    if ~isempty(varargin) %WHistFlag == 1
        if varargin{1} == 1
            WResStruct.WHistCell = WHistCell;
        end
    end
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

function [ val,grad, H ] = Z0_termFunc( Y, Z0, rho, res1, scaleFac, vNum, coordXY )
%Z0_termFunc the terms consisted of z0 in ADMM formulation
eZ0 = exp(Z0);
%     fprintf( 'll = %g, pp = %g\n', full(scaleFac*sum( ( -y.*z0 + eZ0 ))), full(rho / 2 * sum( (z0 +res1).^2 )));
val = sum( scaleFac*( -Y.*Z0 + eZ0 ) ) + rho / 2 * sum( (Z0 +res1).^2 );
grad = scaleFac*(-Y + eZ0) + rho* ( Z0 + res1 );
tmp = scaleFac*eZ0 + rho;
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


function val = LLfunc( Y, W, W0, D, scaleFactor)
    [sLen, ~] = size( Y );
    preY = D(:, :)*W(:, :) + repmat( W0, sLen, 1 );
    val = -scaleFactor*Y(:, :).*preY(:, :) + exp( preY(:, :) );
    val = sum(val(:));
end

