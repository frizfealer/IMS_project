% function theta = drchrnd(alpha,n)
%
% Generate n samples from Dirichlet distribution with 
% vector of parameters alpha (d x 1)
% 
% Note: Each Dirichlet sample is obtained by 
%       getting d independent samples from Gamma and normalizing them to sum to 1
%
% Uses gamrnd and rndcheck functions from the Matlab Statistics Toolbox
%
% josef@robots.ox.ac.uk
% 21/7/2004
function theta = drchrnd(alpha,n)

p = length(alpha);
r = zeros(p,n);

if size(alpha,2)>size(alpha,1)
   alpha = alpha';
end;   

if 0   
   for i = 1:n
      theta(:,i) = gamrnd(alpha,1);
      theta(:,i) = theta(:,i) / sum(theta(:,i));
   end;       
else % faster version
   theta = gamrnd(repmat(alpha,1,n),1,p,n);   
   theta = theta ./ repmat(sum(theta,1),p,1);
end;


return;