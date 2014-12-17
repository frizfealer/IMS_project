% function [Pw_z,Pd_z,Pz,Li] = do_pLSA(alpha,beta,N,nd,Par)
%
% 1. Generate a corpus of documents using the latent Dirichlet allocation (LDA) model
%    with parameters alpha and topics beta.
%
% 2. Fit the pLSA model to the sampled collection
%
% Input:
%        alpha (K x 1) parameters of the Dirichlet distribution
%        beta  (M x K) probability of a word given topic 
%        M     (1 x 1) vocabulary size
%        K     (1 x 1) number of topics
%        N     (1 x 1) number of words per document
%        nd    (1 x 1) number of documents
%
% Output:  
% Li     ... likelihood for each pLSA iteration
% Parameters of the pLSA model:
%   Pz   ... P(z)
%   Pd_z ... P(d|z) 
%   Pw_z ... P(w|z) 
%
% Josef Sivic
% josef@robots.ox.ac.uk
% 30/7/2004
function [Pw_z,Pd_z,Pz,Li] = do_pLSA(alpha,beta,N,nd,Par)
 
% defualt parameters
if nargin < 1 
   N  = 100;  % number of words per document
   nd = 1000; % number of documents   
   alpha = [.2 .2 .2]';   
   beta  = [.20 .20 .20 .20 .20   0   0   0   0   0   0   0;
      0  0   0  .20 .20 .20 .20 .20   0   0   0   0;
      0  0   0   0   0   0   0    0  .25 .25 .25 .25]';

   Par.doplot = 0;
   Par.maxit  = 100;
   Par.Leps   = 1;
end;
        
% vocabulary size, number of topics
M = size(beta,1);
K = size(beta,2);

% generate corpus
X = sample_documents(alpha,beta,N,nd);

% fit the pLSA model
tic;
[Pw_z,Pd_z,Pz,Li] = pLSA_EM(X,K,Par);
fprintf('\n %d Iterations of EM took %.2f sec, Likelihood change:%f  \n',length(Li),toc,Li(end)-Li(end-1));

% plot likelihood of each iteration
if Par.doplot >= 1
   ff(1)=figure(1); clf;
   plot(Li,'b.-');
   %set(ff(1),'Position',[1137 772 452 344]);
   title('Log-likelihood'); xlabel('Iteration'); ylabel('Log-likelihood');
   grid on;
end;

% show results      
if 1      
   fprintf('\n True topics \n');
   beta
   
   fprintf('\n \n Estimated topics \n');
   Pw_z
end;   




