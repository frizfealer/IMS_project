function [ output_args ] = generative_model_MG( uVec, sVec, piVec )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
gridSize= [30, 30];
mu1 = [3 3];
SIGMA1 = [2 0; 0 2];
rng('default');  % For reproducibility
pi = [0.7 0.3];
num = 1000;
mu2 = [8 8];
SIGMA2 = [3 0; 0 3];
can = zeros(100, 1);
n=floor(100*pi(1));
can(1:n) = 1;
for i = 2:length(pi)
    nn=floor(100*pi(i));
    if i == length(pi)
        can(n+1:end)=i;
    else
        can(n+1:(n+1+nn)) = i;
    end
    n = n+1+nn;
end
for i = 1:num
    ins = randi(100);
    if can(ins) == 1
        r3(i,:) = mvnrnd(mu1,SIGMA1,1);
    elseif can(ins) == 2
        r3(i,:) = mvnrnd(mu2,SIGMA2,1);
    end
end

% rng('default');  % For reproducibility
% r2 = mvnrnd(mu2,SIGMA2,100);
plot(r3(:,1),r3(:,2),'+');
[label model L]=vbgm(r3',10); 
end

