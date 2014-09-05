function [ pureAIC, AICC, AICC_od, c_hat ] = computeAICC( LL, LLY, k, n )
%compueAICC compute AIC corrected function
% overDisFlag: over-dispersion flag for Poisson regression
% n: sample number, in our case, is the size(Y, 1)*size(Y, 2)
% k: variable number
% LL: log-likelihood of the input model
% LLY: log-likelihood of the data 
% AICC_od: the ACC with over dispersion modification
pureAIC = -2*LL + 2*k;
AICC = -2*LL + 2*k + ( 2*k*(k+1) ) / max( 1, (n-k-1) );
% logY = log(Y);
% logY(isinf(logY)) = 0;
%LLY = LLFunc( Y, [], [], logY );
resDev = 2*( abs( LLY- LL ) );
c_hat = resDev / max(1, (n-k-1));
AICC_od = -2*LL/c_hat + 2*k + ( 2*k*(k+1) ) / max( 1, (n-k-1) );

end

