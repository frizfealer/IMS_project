function [ val, meanVal, stdVal ] = LP_DL_Poiss( aMatrix, Y, W, W0, D, lambda, phi, theta, scaleFactor, logFY, varargin )
%LP_DL_Poiss Compute the negative log posterior for
%Dictioanry Learning with Poisson distribution
% aMatrix [h w], an indicator matrix of with 1 means for tranining
% Y [s h w], W [m h w] or W[m h*w], W0 [h, w] or W0 [h*w, 1], D [s m]
% lambda, phi, theta, scalar
% logFY [s h w], log of the factorial of Y, can be empty
% scaleFactor [s h w], the scaling factor (0~1) for the log-likelihood
% terms, if empty, then the value is 1
% shoud minimize this value
% meanFlag:, is == 1, compute mean of LL and std of LL for input samples
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
firstTwoTermsMat = -Y(:, idx).*preY(:, idx) + exp( preY(:, idx) );
if ~isempty( logFY )
    firstTwoTermsMat = firstTwoTermsMat - logFY(:, idx);
end
if isempty( scaleFactor )
    scaleFactor = 1;
end
firstTwoTermsMat = firstTwoTermsMat * scaleFactor;
firstTwoTerms = sum( sum( firstTwoTermsMat ) );
% firstTwoTerms sum( Y(:).*z0(:) - exp(z0(:)) ) 
val = firstTwoTerms + lambda * norm( W(:), 1 ) ...
    + phi * norm( D(:), 1 ) + theta * ( norm( Whminus(:), 2 ) + norm( Wwminus(:), 2 ) ); %theta * ( norm( Whminus(:), 1 ) + norm( Wwminus(:), 1 ) );
if length(varargin) == 1 && varargin{1} == 1 % meanFlag == 1
    meanVal = mean(firstTwoTermsMat(:) + lambda * norm( W(:), 1 ) ...
    + phi * norm( D(:), 1 ) + theta * ( norm( Whminus(:), 1 ) + norm( Wwminus(:), 1 ) ) );
    stdVal = std(firstTwoTermsMat(:)) + lambda * norm( W(:), 1 ) ...
    + phi * norm( D(:), 1 ) + theta * ( norm( Whminus(:), 1 ) + norm( Wwminus(:), 1 ) );
else 
    meanVal = [];
    stdVal = [];
end