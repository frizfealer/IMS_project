function [ val ] = validationOnTesting( aMatrix, inY, inW, inW0, inD )
%validationOnTesting using the trained W, W0, and D on testing data
%inY [s h w] inW [m h w] inW0 [h w] inD [s m]
%aMatrix [h w], with 1 on grid means testing, 0 means trainning
[sLen, ~, ~] = size( inY );
preY = inD* inW(:,:) + repmat( inW0(:)', sLen, 1); %preY [s ,w*h]
lambda = floor(exp(preY));
idx = aMatrix == 1;
val = inY(:, idx).*log( lambda(:, idx) ) - lambda(:,idx);
val = sum( sum( val ) );
% for i = 1:hei
%     for j = 1:wid
%         if aMatrix(i,j) == 1 && dataMask2(i,j) == 1
%             val = val + inY(:,i+(j-1)*hei)'*log( lambda(:,i+(j-1)*hei) ) - sum( lambda(:,i+(j-1)*hei) ) - sum(logFInsBY(:,i+(j-1)*hei));
%         end
%     end
% end

end

