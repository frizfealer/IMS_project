function [ val, preY ] = validationOnTesting( aMatrix, Y, D, W, W0, lambda, linkFunc, varargin )
%validationOnTesting using the trained W, W0, and D on testing data Y
%--------------------------------------------------------------------------
% validationOnTesting: validation on the leave-out data
%--------------------------------------------------------------------------
% DESCRIPTION:
%   using trained D on leave-out Y.
%
% INPUT ARGUMENTS:
%   aMatrix, the indicator matrix
%	Y, the input data cube
%	D, the trained dictionary
%	linkFunc, should be 'identity', under construction
% OUTPUT ARGUMENTS:
%   val, log-posterior value,
%	preY, prediciton Y

[sLen, ~, ~] = size( Y );
BlkDS = conBLKDS( Y );
tIdx = BlkDS.indMap==1 & aMatrix == 0;
preY = D*W(:,tIdx) + repmat( W0( tIdx )', sLen, 1 );
testY = Y(:, tIdx);
% if strcmp( linkFunc, 'identity' ) == 1
%     val = testY.*log(preY) - preY;
%     val(isnan(val)) = 0;
%     val(val==-inf) = 0;
%     val = sum(val(:));
% elseif strcmp( linkFunc, 'log' ) == 1
%     val = testY.*preY - exp(preY);
%     val = sum(val(:));
% elseif strcmp( linkFunc, 'log_gaussian' ) == 1
%     ttt = log(testY);
%     ttt(ttt==-inf)=0;
%     val = -(ttt-preY).^2;
%     val = sum(val(:));
% elseif strcmp( linkFunc, 'negative_binomial' ) == 1
%     kappa = varargin{1}; 
%     val = testY*log(kappa)+testY.*preY ...
%         -(testY+1/kappa).*log( 1+kappa*exp(preY) ) ...
%         + gammaln( testY+1/kappa ) - gammaln( 1/kappa );
%     val = sum(val(:));
% end



% %compute scaleFactor
% scaleFactor =  1e-2;
ins = max( D, [], 1);
mLen = length( find( ins >= 0.01 ) );
if ~isempty(mLen)
    [WResStruct] = updateW_ADMM_testing( testY, D, 200, [], [], 1, lambda, 5e-3, 1e-2, [] );
    newWInfo = [];
%     for i = 1:size(testY(:,:), 2)
%         newWInfo{i} = find( WResStruct.W(:, i) > max(WResStruct.W(:, i))*1e-2);
%     end
%     [WResStruct] = updateW_ADMM_testing( testY, D, 200, [], [], scaleFactor, lambda, 5e-3, 1e-2, newWInfo );
    preY = D* WResStruct.W(:, :) + repmat( WResStruct.W0(:)', sLen, 1 ); 
    idx = find(preY ~= 0);
    val = testY(idx).*log(preY(idx)) - preY(idx);
    val = sum(val(:));
else
    val = 0;
end
end

