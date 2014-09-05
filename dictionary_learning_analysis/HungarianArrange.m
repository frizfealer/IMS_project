function [ mfitD, Matching ] = HungarianArrange( fitD, gD )
%HungarianArrange Use Hungarian match to rearrange fitD.
for i = 1:size(gD, 2)
    gD(:, i) = gD(:,i) / norm( gD(:,i), 2 );
    fitD(:, i) = fitD(:, i ) / norm( fitD(:, i), 2 );
end
costMat1 = gD' * fitD;
%% after Hungarian matching
[Matching, cost]=Hungarian( -costMat1 );
mfitD = zeros( size( fitD, 1 ), size( fitD, 2 ) );
%modify pD to original matrix by Hungarian results
for i=1:size( Matching, 1 )
    mfitD(:,i) = fitD(:,find(Matching(i,:)==1));
end


end

