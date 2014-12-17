% function bin = multinom_rnd(theta)
%
% Sample from a histogram (multinomial distribution) 
% with bin probabilities theta, returns sampled bin number
%
% josef@robots.ox.ac.uk
% 21/7/2004
function bin = multinom_rnd(theta)

cmtheta = cumsum(theta);
theta_ranges = [[0; cmtheta(1:end-1)] [cmtheta(1:end)]];
rn = rand(1);
bin = find(rn>theta_ranges(:,1) & rn<=theta_ranges(:,2));

return;
