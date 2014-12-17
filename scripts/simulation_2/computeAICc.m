function [ AICC_3 ] = computeAICc( gY, LL, D_pLSA, BlkDS, k)
%W_pLSA = D_pLSA \ gY(:, BlkDS.indMap == 1);
samNum = length( gY(1, BlkDS.indMap == 1 ) ) * size(gY, 1);

%LL = -( sum(sum(gY(:,BlkDS.indMap==1)-D_pLSA * W_pLSA)) );
m =  k*( size(gY, 1) + length( gY(1, BlkDS.indMap == 1 ) ) );
n = samNum;
AICC_3 = -2*LL + 2*m + ( 2*m*(m+1) ) / max( 1, (n-m-1) );


end

