function [ val, meanVal, stdVal ] = LP_DL_Poiss( LINK_FUNC, aMatrix, Y, W, W0, D, lambda, phi, theta, scaleFactor, logFY, varargin )
%--------------------------------------------------------------------------
%LP_DL_Poiss: ompute the negative log posterior for Dictioanry Learning 
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
%   meanVal, stdVal, mean and standard deviation of log posterior of the
%   input.

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
% firstTwoTerms = 0;
% for i = 1:hei
%     for j = 1:wid
%         if aMatrix(i,j) == 1
%             firstTwoTerms = firstTwoTerms - ( scaleFactor(:, j+(i-1)*hei).*Y(:,i+(j-1)*hei) )'*preY(:,i+(j-1)*hei) + ...
%                 sum( scaleFactor(:, j+(i-1)*hei).*exp(preY(:,i+(j-1)*hei)) );
%         end
%     end
% end
idx = aMatrix == 1;
if strcmp( LINK_FUNC, 'log' ) == 1
    firstTwoTermsMat = -Y(:, idx).*preY(:, idx) + exp( preY(:, idx) );
    if ~isempty( logFY )
        firstTwoTermsMat = firstTwoTermsMat - logFY(:, idx);
    end
elseif strcmp( LINK_FUNC, 'identity' ) == 1
    %     preY(preY==0)=1e-32;
    tmp = log( preY(:, idx) + 1e-8 );
    firstTwoTermsMat = -Y(:, idx).*tmp + preY(:, idx);
    if ~isempty( logFY )
        firstTwoTermsMat = firstTwoTermsMat - logFY(:, idx);
    end
elseif strcmp( LINK_FUNC, 'log_gaussain' ) == 1
    firstTwoTermsMat = (log(Y(:, idx)+1e-32)-preY(:, idx)).^2;
elseif strcmp( LINK_FUNC, 'negative_binomial' ) == 1
    kappa = 1e-2;
    if length(varargin) == 2 && ~isempty( varargin{2} )%set varargin{1} as []
        kappa = varargin{2};
    end
    dataTerms = -Y(:, idx)*log(kappa) - Y(:, idx).*preY(:, idx) ...
        + (Y(:, idx)+1/kappa).*log( 1+kappa*exp(preY(:, idx)) ) ...
        - gammaln( Y(:, idx)+1/kappa ) + gammaln( 1/kappa );
    firstTwoTermsMat = dataTerms;
end
if isempty( scaleFactor )
    scaleFactor = 1;
end
firstTwoTermsMat = firstTwoTermsMat * scaleFactor;
firstTwoTerms = sum( sum( firstTwoTermsMat ) );
% firstTwoTerms sum( Y(:).*z0(:) - exp(z0(:)) ) 
val = firstTwoTerms + lambda * norm( W(:), 1 ) ...
    + phi * norm( D(:), 1 ) + theta * ( norm( Whminus(:) )^2 + norm( Wwminus(:) )^2 ); %theta * ( norm( Whminus(:), 1 ) + norm( Wwminus(:), 1 ) );
if length(varargin) == 1 && varargin{1} == 1 % meanFlag == 1
    meanVal = mean(firstTwoTermsMat(:) + lambda * norm( W(:), 1 ) ...
    + phi * norm( D(:), 1 ) + theta * ( norm( Whminus(:), 1 ) + norm( Wwminus(:), 1 ) ) );
    stdVal = std(firstTwoTermsMat(:)) + lambda * norm( W(:), 1 ) ...
    + phi * norm( D(:), 1 ) + theta * ( norm( Whminus(:), 1 ) + norm( Wwminus(:), 1 ) );
else 
    meanVal = [];
    stdVal = [];
end