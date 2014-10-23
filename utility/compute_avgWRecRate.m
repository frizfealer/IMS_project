function [ val ] = compute_avgWRecRate( W1, W2 )
%compute_avgWRecRate compute average W recovery rate according to my paper
%W1 is the ground truth abundance, where as W2 is the learned abundance
val = mean(1-abs(W1-W2)./W1);


end

