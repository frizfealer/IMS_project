function [D, W, molD, insW] = postProcessResults( D, W, thresD, thresW )
%postProcessResults remove D and W that has smaller signals then thresD and
%thresW
%thresD: percentage of the largest value in D, default == 0.01
%thresW: percentage of the largest value in W in each grid cell, default shoud be 0.01
if isempty( thresD )
    thresD = 0.01;
end
if isempty( thresW )
    thresW = 0.01;
end
%% on D
[~, mLen] = size(D);
%remove small signals in each dictionary elements
for i = 1:mLen
    cD = D(:, i);
    maxD = max(cD);
    D(cD<maxD*thresD, i) = 0;
end
%remove dictionary elements with small signals.
ins = zeros(mLen, 1);
for i = 1:mLen
ins(i) = max(D(:,i));
end
molD = find( ins>1*thresD );
% D = D(:, molD);

%% on W
for i = 1:size(W, 2)
    maxW = max(W(:, i));
    cW = W(:,i);
    W(cW < maxW*thresW, i) = 0;
end

molI = molD;
W = W(molI, :);
insW = zeros( length(molI), 1 );
for i = 1:length(insW)
    insW(i) = max(W(i,:));
end
% molW = find( ins>ins*thresW );
% molI = intersect( molW, molD );
molI = molD;
D = D(:, molI);
end

