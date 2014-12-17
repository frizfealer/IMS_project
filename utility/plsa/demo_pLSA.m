% Probabilistic latent semantic analysis (pLSA) demo
%
%   1. Sample a corpus of documents from the latent Dirichlet allocation model (LDA) [1]
%   2. Fit the probabilistic latent semantic analysis model (pLSA)  [2,3]
%
% This demo code is based on the paper:
%  J. Sivic, B. C. Russell, A. Efros, A. Zisserman and W. T. Freeman,
%  Discovering objects and their location in images, ICCV 2005.
%  http://www.robots.ox.ac.uk/~vgg/publications/papers/sivic05b.pdf 
%
% References:
% [1] D. Blei, A. Ng, and M. Jordan: Latent Dirichlet allocation, Journal
% of Machine Learning Research, 3:993-1022, 2003.
% [2] T. Hofmann: Probabilistic Latent Semantic Analysis, 
% Proc. of the 15th Conf. on Uncertainty in Artificial Intelligence (UAI'99) 
% [3] T. Hofmann: Unsupervised Learning by Probabilistic Latent Semantic
% Analysis, Machine Learning Journal, 42(1), 2001, pp.177.196 
%
%
% The LDA corpus sampling code is based on code by
% Bryan Russell <brussell@csail.mit.edu>
%
% Josef Sivic 
% josef@robots.ox.ac.uk
% 30/7/2004

% Set plSA EM parameters
Par.maxit  = 200; % maximum number of iterations of EM
Par.Leps   = 0.1; % EM stopping condition: minimum change in log-likelihood
Par.doplot = 1;   % Ploting verbosity

% Corpus parameters
N     = 10;   % number of words per document
nd    = 1000; % number of documents
alpha = [.2 .2 .2]'; % Dirichlet parameter

% Topic distributions: (3 topics, vocabulary of 12 distinct terms)
beta  = [.25 .25 .25 .25   0   0   0   0   0   0   0   0;
           0  0   0   0  .25 .25 .25 .25   0   0   0   0;
           0  0   0   0   0   0   0   0  .25 .25 .25 .25]';

% run pLSA demo       
[Pw_z,Pd_z,Pz,Li] = do_pLSA(alpha,beta,N,nd,Par);

