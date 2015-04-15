% function [Pw_z,Pd_z,Pz,Li] = pLSA_EM(X,Fixed_Pw_z,K,Learn)
%
% Probabilistic Latent semantic alnalysis (pLSA)
%
% Notation:
% X ... (m x nd) term-document matrix (observed data)
%       X(i,j) stores number of occurrences of word i in document j
%
% m  ... number of words (vocabulary size)
% nd ... number of documents
% K  ... number of topics
%
% Fixed_Pw_z ... fixed Pw_z density (for recognition only)
%                leave empty in learning
%                N.B. Must be of size m by K
%
% Learn ... Structure holding all settings 
%
% Li   ... likelihood for each iteration
% Pz   ... P(z)
% Pd_z ... P(d|z) 
% Pw_z ... P(w|z) corresponds to beta parameter in LDA
%
% Pz_wd ... P(z|w,d) posterior on z
%
% 
% References: 
% [1] Thomas Hofmann: Probabilistic Latent Semantic Analysis, 
% Proc. of the 15th Conf. on Uncertainty in Artificial Intelligence (UAI'99) 
% [2] Thomas Hofmann: Unsupervised Learning by Probabilistic Latent Semantic
% Analysis, Machine Learning Journal, 42(1), 2001, pp.177.196 
%
% Josef Sivic
% josef@robots.ox.ac.uk
% 30/7/2004
%
% Extended by Rob Fergus
% fergus@csail.mit.edu
% 03/10/05

function [Pw_z,Pd_z,Pz,Li, fL] = pLSA_EM_sparse(X,K, p1, p2, Pz_dw_init)

%% small offset to avoid numerical problems  
ZERO_OFFSET = 1e-7;

%%% Default settings
   Learn.Max_Iterations  = 1000;
   Learn.Min_Likelihood_Change   = 1;   
   Learn.Verbosity = 0;

if Learn.Verbosity
   figure(1); clf;
   title('Log-likelihood');
   xlabel('Iteration'); ylabel('Log-likelihood');
end;

m  = size(X,1); % vocabulary size
nd = size(X,2); % # of documents

% initialize Pz, Pd_z,Pw_z
[~,Pd_z,Pw_z] = pLSA_init(m,nd,K);
Pz_dw = Pz_dw_init;
%% using the original pLSA to estimate C and D
C = zeros(m, K);
D = zeros(K, nd);
alphak = zeros(K, 1);
betai = zeros(K, 1);
for k = 1:K
    C(:,k) = sum(X .* Pz_dw(:,:,k),2);
    alphak(k) = p1*sum( C(:, k) );
    D(k, :) = sum(X.* Pz_dw(:,:,k),1)';
end; 
for i = 1:nd
    betai(i) = p2*sum( D(:, i) );
end
%% run the sprase pLSA
Pz   = ones(K,1)/K; % uniform prior on topics
maxit = Learn.Max_Iterations;

Li = zeros(maxit, 1);
% EM algorithm
for it = 1:maxit   
   fprintf('Iteration %d ',it);
   
   % E-step
   Pz_dw = pLSA_Estep(Pw_z,Pd_z,Pz);
   
   % M-step
   [Pw_z,Pd_z] = pLSA_Mstep_sparse( X, Pz_dw, alphak, betai );

    
   % Evaluate data log-likelihood
   Li(it) = pLSA_logL(X,Pw_z,Pz,Pd_z, alphak, betai);   
        
   % plot loglikelihood
   if Learn.Verbosity>=3
      figure(ff(1));
      plot(Li,'b.-');
   end;
      
   %%% avoid numerical problems.
   Pw_z = Pw_z + ZERO_OFFSET;
   
   % convergence?
   dLi = 0;
   if it > 1
     dLi    =  Li(it) - Li(it-1) ;
     if abs(dLi) < Learn.Min_Likelihood_Change,
         fL = sum(sum(X .* log(Pw_z * diag(Pz) * Pd_z' + 1e-8)));
         break; 
     end   
   end;
   fprintf('dLi=%f \n',dLi);
end;
fL = sum(sum(X .* log(Pw_z * diag(Pz) * Pd_z' + 1e-8)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize conditional probabilities for EM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pz,Pd_z,Pw_z] = pLSA_init(m,nd,K)
% m  ... number of words (vocabulary size)
% nd ... number of documents
% K  ... number of topics
%
% Pz   ... P(z)
% Pd_z ... P(d|z)
% Pw_z ... P(w|z)

Pz   = ones(K,1)/K; % uniform prior on topics

% random assignment
Pd_z = rand(nd,K);   % word probabilities conditioned on topic
C    = 1./sum(Pd_z,1);  % normalize to sum to 1
Pd_z = Pd_z * diag(C);

% random assignment
Pw_z = rand(m,K);
C    = 1./sum(Pw_z,1);    % normalize to sum to 1
Pw_z = Pw_z * diag(C);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) E step compute posterior on z,  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pz_dw = pLSA_Estep(Pw_z,Pd_z,Pz)
[m, K] = size(Pw_z);
[nd] = size(Pd_z, 1);
Pz_dw = zeros(m,nd,K);
for k = 1:K
    Pz_dw(:,:,k) = Pw_z(:,k) * Pd_z(:,k)';% * Pz(k);
end;
C = sum(Pz_dw,3);

% normalize posterior
for k = 1:K
    Pz_dw(:,:,k) = Pz_dw(:,:,k) .* (1./C);
end;

end

function [Pw_z,Pd_z] = pLSA_Mstep_sparse(X, Pz_dw, alphak, betai)
[M, ND, K] = size(Pz_dw);
C = zeros(M, K);
D = zeros(K, ND);
Pw_z = zeros( size( C ) );
Pk_d = zeros( size( D ) );
eps = 1e-8;

for k = 1:K
    C(:,k) = sum(X .* Pz_dw(:,:,k),2);
    Pw_z(:, k) = ( C(:, k) - alphak(k) ) ./ ( sum( C(:, k) ) - M*alphak(k) );
    negIdx = find( Pw_z(:, k) < eps );
    posIdx = find( Pw_z(:, k) >= eps );
    if ~isempty( negIdx )
        Pw_z(negIdx, k) = eps;
        n = M - length(negIdx);
        Pw_z(posIdx, k) = ( C(posIdx, k) - alphak(k) ) ./ ( sum( C(posIdx, k) ) - n*alphak(k) );
    end
end

for k = 1:K
    D(k, :) = sum(X.* Pz_dw(:,:, k),1);
end
for i = 1:ND
    Pk_d(:, i) = ( D(:, i) - betai(i) ) ./ ( sum( D(:, i) ) - K*betai(i) );
    negIdx = find( Pk_d(:, i) < eps );
    posIdx = find( Pk_d(:, i) >= eps );
    if ~isempty( negIdx )
        Pk_d(negIdx, i) = eps;
        n = K - length(negIdx);
        Pk_d(posIdx, i) = ( D(posIdx, i) - betai(i) ) ./ ( sum( D(posIdx, i) ) - n*betai(i) );
    end
end;
Pd_z = Pk_d';

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data log-likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = pLSA_logL(X,Pw_z,Pz,Pd_z, alphak, betai)
L = sum(sum(X .* log(Pw_z * diag(Pz) * Pd_z' + 1e-8)));
L = L - sum( log(Pw_z), 1 )*alphak - sum( log(Pd_z), 2)'*betai;
end 





