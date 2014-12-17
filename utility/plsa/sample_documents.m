% function X = sample_documents(alpha,beta,nw,nd)
% 
% Sample documents from the latent Dirichlet allocation model
%
% Input: alpha (k x 1) parameters of the Dirichlet distribution
%        beta  (m x k) probability of a word given topic (m is vocabulary size)
%        nw            number of words per document
%        nd            number of documents
% 
% Output: term document matrix X (m x nd)
%
% josef@robots.ox.ac.uk
% 28/7/2004

function X = sample_documents(alpha,beta,nw,nd)

fprintf('Generating  %d documents \n',nd);
% sample all documents in the corpus

% sample topic weights for each document
theta_v = drchrnd(alpha,nd);

m = size(beta,1); % number of distinct terms
D = zeros(m,nd);
for j = 1:nd

   % sample topic weights for each document
   theta = theta_v(:,j);
   
   % sample all words in a document
   for i = 1:nw
      % sample topic from multinomial with weights theta
      zi = multinom_rnd(theta);
      
      % sample word from the topic
      wi = multinom_rnd(beta(:,zi));
      
      D(wi,j) = D(wi,j)+1;
   end;
   
end;
fprintf('Done \n',nd);
X = D;

return;