function [ val, dataTermsVal, lPenaltyVal, pPenaltyVal, tPenaltyVal ] = LP_DL_Poiss( LINK_FUNC, aMatrix, Y, W, W0, D, lambda, phi, theta, scaleFactor, logFY, varargin )
%--------------------------------------------------------------------------
%LP_DL_Poiss: compute the negative log posterior for Dictioanry Learning 
%with the Poisson distribution
%--------------------------------------------------------------------------
% DESCRIPTION:
%   In current version, this function only computes the negative log
%   posterior under the Poisson distribution setting. This value should be
%   minimized with dictionary learning.
%
% INPUT ARGUMENTS:
%   aMatrix [h w], an indicator matrix of with 1 means for tranining
%   Y [s h w], W [m h w] or W[m h*w], W0 [h, w] or W0 [h*w, 1], D [s m]
%   lambda, phi, theta, scalar
%   logFY [s h w], log of the factorial of Y, can be empty
%   scaleFactor [s h w], the scaling factor (0~1) for the log-likelihood
%   terms, if empty, then the value is 1.
%   meanFlag, additional variables, if it is == 1, 
%   compute mean of LL and std of LL for input samples.
% OUTPUT ARGUMENTS:
%   val, total negative log posterior in the model.
%   dataTermsVal, lPenaltyVal, pPenaltyVal, tPenaltyVal: the values of each
%   terms

[sLen, hei, wid] = size( Y );

preY = D* W(:, :) + repmat( W0(:)', sLen, 1 ); %preY [s ,w*h]
%compute w(i,j) - w(i-1,j)
%Whminus, a matrix records W_{i,j} - W{i-1,j}
Whminus = zeros( size(W) );
for i = 1:wid
    Whminus(:, 1+(i-1)*hei) = 0;
end
for i = 2:hei
    for j = 1:wid
        Whminus(:, i+(j-1)*hei) = W(:, i+(j-1)*hei) - W(:,i-1+(j-1)*hei);
    end
end
%Wwminus, a matrix records W_{i,j} - W{i,j-1}
Wwminus = zeros( size(W) );
for i = 1:hei
    Wwminus(:,i) = 0;
end
for i = 2:wid
    for j = 1:hei
        Wwminus(:, j+(i-1)*hei) = W(:,j+(i-1)*hei) - W(:,j+(i-2)*hei);
    end
end

%compute the Log posterior
idx = aMatrix == 1;
if strcmp( LINK_FUNC, 'log' ) == 1
    dataLikeliTerms = -Y(:, idx).*preY(:, idx) + exp( preY(:, idx) );
    if ~isempty( logFY )
        dataLikeliTerms = dataLikeliTerms + logFY(:, idx);
    end
elseif strcmp( LINK_FUNC, 'identity' ) == 1
    tmp = preY(:, idx);
    tmp(tmp~=0) = log( tmp(tmp~=0) );
    dataLikeliTerms = -Y(:, idx).*tmp + preY(:, idx);
    if ~isempty( logFY )
        dataLikeliTerms = dataLikeliTerms + logFY(:, idx);
    end
elseif strcmp( LINK_FUNC, 'log_gaussain' ) == 1
    dataLikeliTerms = (log(Y(:, idx))-preY(:, idx)).^2;
elseif strcmp( LINK_FUNC, 'negative_binomial' ) == 1
    kappa = 1e-2;
    if length(varargin) == 2 && ~isempty( varargin{2} )%set varargin{1} as []
        kappa = varargin{2};
    end
    dataTerms = -Y(:, idx)*log(kappa) - Y(:, idx).*preY(:, idx) ...
        + (Y(:, idx)+1/kappa).*log( 1+kappa*exp(preY(:, idx)) ) ...
        - gammaln( Y(:, idx)+1/kappa ) + gammaln( 1/kappa );
    dataLikeliTerms = dataTerms;
end
if isempty( scaleFactor )
    scaleFactor = 1;
end
dataLikeliTerms = dataLikeliTerms * scaleFactor;
dataLikeliVal = sum( sum( dataLikeliTerms ) );
% firstTwoTerms sum( Y(:).*z0(:) - exp(z0(:)) ) 
val = dataLikeliVal + lambda * norm( W(:), 1 ) ...
    + phi * norm( D(:), 1 ) + theta * ( norm( Whminus(:) )^2 + norm( Wwminus(:) )^2 ); %theta * ( norm( Whminus(:), 1 ) + norm( Wwminus(:), 1 ) );
dataTermsVal = dataLikeliVal;
lPenaltyVal = lambda * norm( W(:), 1 ); 
pPenaltyVal = phi * norm( D(:), 1 );
tPenaltyVal =  theta * ( norm( Whminus(:) )^2 + norm( Wwminus(:) )^2 );
end