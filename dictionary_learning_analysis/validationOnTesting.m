function [ val, preY ] = validationOnTesting( aMatrix, Y, D, W, W0, lambda, linkFunc )
%validationOnTesting using the trained W, W0, and D on testing data Y

[sLen, ~, ~] = size( Y );
BlkDS = conBLKDS( Y );
tIdx = BlkDS.indMap==1 & aMatrix == 0;
preY = D*W(idx,tIdx) + repmat( W0( tIdx )', sLen, 1 );
testY = Y(:, tIdx);
if strcmp( linkFunc, 'identity' ) == 1
    val = testY.*log(preY+1e-32) - preY;
    val = sum(val(:));
elseif strcmp( linkFunc, 'log' ) == 1
    val = testY.*preY - exp(preY);
    val = sum(val(:));
end



% %compute scaleFactor
% scaleFactor =  1e-2;
% if mLen~=0
%     [WResStruct] = updateW_ADMM_testing( testY, D, 200, [], [], scaleFactor, lambda, 5e-3, 1e-2, [] );
%     newWInfo = [];
%     for i = 1:size(testY(:,:), 2)
%         newWInfo{i} = find( WResStruct.W(:, i) > max(WResStruct.W(:, i))*1e-2);
%     end
%     [WResStruct] = updateW_ADMM_testing( testY, D, 200, [], [], scaleFactor, lambda, 5e-3, 1e-2, newWInfo );
%     preY = D* WResStruct.W(:, :) + repmat( WResStruct.W0(:)', sLen, 1 ); 
%     val = testY.*preY - exp( preY);
%     val = sum(val(:));
% else
%     val = 0;
% end
end

