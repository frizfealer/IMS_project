function [WResStruct] = updateW_ADMM_testing( Y, D, itNum, logFY, initVar, scaleFactor, lambda, varargin )
%--------------------------------------------------------------------------
% updateW_ADMM_testing: update W with ADMM formulation:
% z0 = D*W, z1 = W
%--------------------------------------------------------------------------
% DESCRIPTION:
%   This function is for testing the trained dictionary on the leave-out data.
%	This function only has a hyperparameter: lambda.
%
% INPUT ARGUMENTS:
%   Y, size of [s n], s is # m/z and n is # samples.
%   D, the trained dictionary
%   itNum, iteration number
%   logFY, log factorial Y, can be []
%	initVar, initial values of W, can be []
%	scaleFactor, should be 1
%	lambda, hyperparameter
%	Wtol, additional variable 1, the tolerance of W difference 
%	default is 5e-3
%	D_LOWER_BOUND, additional variable 2, the smallest value of D, default is 1e-2
%	newWInfo, additional variable 3, under construction
% OUTPUT ARGUMENTS:
%   WResStruct, a data structre includes W, W0, z0, z1

%% set variables
newWInfo = [];
if ~isempty(varargin)  
    if length(varargin) == 2
        Wtol = varargin{1};
        D_LOWER_BOUND = varargin{2};
    else
        Wtol = 5e-3;
        D_LOWER_BOUND = 1e-2;
    end
    if length(varargin) == 3
        Wtol = varargin{1};
        D_LOWER_BOUND = varargin{2};
        newWInfo = varargin{3};
    end
else
    Wtol = 5e-3;
    D_LOWER_BOUND = 1e-2;
end


fprintf( 'iteration number for W_update_ADMM = %d\n', itNum );
[sLen, mLen] = size(D);
[~, nLen] = size( Y(:, :) );
%% initialize variables
if isempty( initVar )
    z0 = log(Y); z0(z0==-inf)=0; %initialize zo as log(Y)
    z0(z0==0)=1e-32;
    z1 = zeros( mLen, nLen );
    W = zeros( mLen, nLen ); W0 = zeros( 1, nLen );
else
    z0 = initVar.z0(:, :);
    z0(z0==0)=1e-32;
    z1 = initVar.z1(:, :);
    W = initVar.W(:, :); W0 = initVar.W0(:, :);
end

LPAryADMM = zeros( itNum+1, 1 );
resRecAry = zeros( itNum, 4, 1 );
curRhoAry = zeros( itNum, 1 );
curRhoAry(1, :) = 1;
u0 = sparse( zeros( size( z0(:, :) ) ) );
u1 = sparse( zeros( size( z1(:, :) ) ) );

moutD = [D ones(sLen, 1)];
tmpEye = eye(mLen+1, mLen+1);
tmpEye(end, end) = 0;
aPartofCforUW = sparse( [moutD; tmpEye] );
CforUW = replicateC( aPartofCforUW, nLen );
if ~isempty(newWInfo)
    [CforUW] =  modifyCWithSparsity( CforUW, newWInfo, 1:nLen, mLen+1, size(CforUW, 1) );
end


%% main optimization process
    LPAryADMM(1) = LLfunc( Y, W, W0, D, scaleFactor);
    for itNumADMM = 1:itNum
%         tic;
        preW = W;
        fprintf( 'LP: %g ', LPAryADMM(itNumADMM) );
        curRho = curRhoAry(itNumADMM);
        %% update W
        fprintf( 'updating W... ' );
        tmpD = [z0+1/curRho*u0; z1+1/curRho*u1; sparse(1, nLen) ];
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
        %% update Z1
        fprintf( 'updating z1...' );
        preZ1 = z1;
        dW = W(:);
        U1 = u1(:);
        C = speye(length(dW));
        d = dW-1/curRho*U1;
        
        opts.init=1;
        opts.x0 = z1(:);
        opts.tFlag=5;       % run  to Xtol
        opts.tol=1e-4;
        opts.nFlag=0;       % without normalization
        opts.rFlag=0;       % the input parameter 'rho' is a ratio in (0, 1)
        opts.lFalg = 1;
        lambda = lambda/curRho;
        opts.rsL2 = 0;
        [tmp, funVal]=LeastR(C, d, lambda, opts);
        z1 = reshape( tmp,  size(z1) );





        %% update u0, u1
        fprintf( 'updating u0... ' );
        preY = D*W(:,:) + repmat( W0(:)', sLen, 1);
        %primal residual r0
        r0 = z0 - preY;
        u0 = u0 + curRho*r0;
        %primal residual r1
        r1 = z1 - W;
        u1 = u1 + curRho*r1;
        %compute the primal residual all
        rf = norm( r0(:) ) + norm( r1(:) );
        pDim = length( r0(:) ) + length( r1(:) );
        %compute the dual variable to see if it is converge
        %dual residual s0
        s0 = curRho*( D'*( z0 - preZ0 ) );
        s1 = curRho*( z1 - preZ1 );
        sf = norm( s0(:) ) + norm( s1(:) );
        nDim = length( s0(:) ) + length( s1(:) );
        
        EPS_ABS = 1e-4;
        EPS_REL = 1e-3;
        
        relEpsPri = max( [norm( z0(:) ) norm( preY(:) ) ...
            norm( z1(:) ), norm( W(:) ) ] );
        epsPri = sqrt(pDim)*EPS_ABS + relEpsPri*EPS_REL;
        
        
        tmp = D'*u0;
        tmp = [tmp; sum(u0)];
        if isempty( tmp )
            tmp = 0;
        end
        maxRel = max( [ norm( tmp(:) ), norm( u1(:) ) ] );        
        epsDual = sqrt(nDim)*EPS_ABS + maxRel*EPS_REL;
        
        resRecAry(itNumADMM, :) = [rf, epsPri, sf, epsDual];
        tmp = full( max( abs( W(:) -  preW(:) ) ) );
        fprintf( '%d: %g %g %g %g %g %g \n',itNumADMM, rf, epsPri, sf, epsDual, curRhoAry(itNumADMM), full(tmp) );
%         time = toc;
%         fprintf( 'time: %g\n', time );        
        %% breaking condition
        if ( rf < epsPri && sf < epsDual ) && ( tmp < Wtol )
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
    WResStruct.z1 = sparse( z1(:, :) );
    WResStruct.u1 = sparse( u1(:, :) );
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
%     fprintf( 'll = %g, pp = %g\n', full(scaleFac*sum( ( -y.*z0 + eZ0 ))), full(rho / 2 * sum( (z0 +res1).^2 )));
val = sum( scaleFac*( -Y.*log(Z0+1e-32) + Z0 ) ) + rho / 2 * sum( (Z0 +res1).^2 );
grad = scaleFac*(-Y./(Z0+1e-32) + 1) + rho* ( Z0 + res1 );
tmp = scaleFac*(Y./(Z0+1e-32).^2) + rho;
H = sparse( coordXY, coordXY, tmp, vNum, vNum );
end

function val = LLfunc( Y, W, W0, D, scaleFactor)
    [sLen, ~] = size( Y );
    preY = D(:, :)*W(:, :) + repmat( W0, sLen, 1 );
    val = -scaleFactor*(Y(:, :).*log(preY(:, :)+1e-32))+ scaleFactor*preY(:, :);
    val = sum(val(:));
end


function [C] = modifyCWithSparsity( CforUW, newWInfo, loc, mLen, nrPart1 )
    for i = 1:length(loc)
        cLoc = loc(i);
        cW = newWInfo{cLoc};
        cW = cW+(i-1)*mLen;
        cWRange = (1+(i-1)*mLen):(i*mLen);
        target = setdiff(cWRange, cW);
        CforUW(1:nrPart1, target) = 0;
    end
    C = CforUW;
end
