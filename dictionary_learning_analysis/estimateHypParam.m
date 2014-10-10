function [ lambda, theta, phi ] = estimateHypParam( eW, eD, DTemplate, BlkDS )
%estimateHypParam estimate hyper parameters
[mLen, hei, wid] = size( eW );
%compute lambda
numSamples = length(find(BlkDS.indMap==1));
lambda = mLen*numSamples / sum( eW(:) );
%compute theta
% nRow = (hei-1)*wid*mLen;
% nCol = (wid-1)*hei*mLen;
Whminus = zeros( size(eW) );
for i = 2:hei
    for j = 1:wid
        Whminus(:, i+(j-1)*hei) = eW(:, i+(j-1)*hei) - eW(:,i-1+(j-1)*hei);
    end
end
Wwminus = zeros( size(eW) );
for i = 2:wid
    for j = 1:hei
        Wwminus(:, j+(i-1)*hei) = W(:,j+(i-1)*hei) - W(:,j+(i-2)*hei);
    end
end
numSamDiff = 0;
[Rall] = genSparseGroupingMatrix( hei, wid, 1 );
for i = 1:BlkDS.blkNum
    Rb = Rall(:, BlkDS.B2GMap{i});
    ins = sum(abs(Rb), 2);
    Rb(ins<2,:) = [];
    numSamDiff = numSamDiff + mLen * size( Rb, 1 );
end
% theta = ( numSamDiff ) / ( sum( abs( Whminus(:) ) ) + sum( abs( Wwminus(:) ) ) );
theta = ( numSamDiff ) / ( norm( Whminus(:), 2 ) + norm( Wwminus(:), 2 ) );
%compute phi
nNonZ = length( find( DTemplate == 1 ) );
phi = nNonZ / sum( eD(:) );

end

