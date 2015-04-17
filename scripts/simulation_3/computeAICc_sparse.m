function [ AICC_3 ] = computeAICc_sparse( samNum, LL, D)
%W_pLSA = D_pLSA \ gY(:, BlkDS.indMap == 1);
%LL = -( sum(sum(gY(:,BlkDS.indMap==1)-D_pLSA * W_pLSA)) );
m =  length( find(D(:)>=1e-2) );
n = samNum;
AICC_3 = -2*LL + 2*m + ( 2*m*(m+1) ) / max( 1, (n-m-1) );


end

