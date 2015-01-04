function [ val ] = validationOnTesting( aMatrix, Y, D )
%validationOnTesting using the trained W, W0, and D on testing data
%inY [s h w] inW [m h w] inW0 [h w] inD [s m]
%aMatrix [h w], with 1 on grid means testing, 0 means trainning

%% preprocessing D
[~, mLen] = size(D);
ins = zeros(mLen, 1);
for i = 1:mLen
ins(i) = max(D(:,i));
end
idx = ins>1e-2;
D = D(:, idx);
[~, mLen] = size(D);

[sLen, ~, ~] = size( Y );
BlkDS = conBLKDS( Y );
tIdx = BlkDS.indMap==1 & aMatrix == 0;

testY = Y(:, tIdx);
testW = zeros( mLen, length(find(tIdx)==1) );
testW0 = zeros( 1, length(find(tIdx)==1) );
%compute scaleFactor
tmp = log(testY); tmp(tmp==-inf)=0; 
scaleFactor =  1 / ( max( testY(:) ) / max( tmp(:) ) );

[WResStruct] = updateW_ADMM_testing( testY, D, 200, [], [], scaleFactor );


preY = D* WResStruct.W(:, :) + repmat( WResStruct.W0(:)', sLen, 1 ); 
val = testY.*preY - exp( preY);
val = sum(val(:));

end

