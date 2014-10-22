function [ AICC_3 ] = computeAICc( gY, D_pLSA, BlkDS, m)
W_pLSA = D_pLSA \ gY(:, BlkDS.indMap == 1);
samNum = length( gY(1, BlkDS.indMap == 1 ) );
LL = -( sum(sum(gY(:,BlkDS.indMap==1)-D_pLSA * W_pLSA)) );
k =  m*samNum+m*samNum;
n = samNum;
AICC_3 = -2*LL + 2*k + ( 2*k*(k+1) ) / max( 1, (n-k-1) );


end

