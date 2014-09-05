function [ pVal ] = computePVal_PoissonReg( LL_p, LL_q, p, q )
%computePVal_PoissonReg Summary of this function goes here
% IGNORE_INT = 1e-2;
% if isempty(p)
%     p = length( find( W1(:) > IGNORE_INT ) );
% end
% if isempty(q)
%     q = length( find( W2(:) > IGNORE_INT ) );
% end
% LL1 = Poisson_LL_Func( Y, D1, W1, [] );
% LL2 = Poisson_LL_Func( Y, D2, W2, [] );
dev = 2*(abs(LL_p-LL_q));
pVal = 1-chi2cdf( dev, abs(p-q) );
end

