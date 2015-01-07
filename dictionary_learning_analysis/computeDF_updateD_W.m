function [ uDIonName, uDTemplate, usedEle, uSpeciesM, uD, uWInfo, dfW ] = computeDF_updateD_W( expRec, DTemplate, DIonName, SpeciesM )
%computeDF_updateD_W compute d.f. of D and W, update D, W
D = expRec.D;
D(D<1e-2) = 0;
uDTemplate = zeros( size(D) );
usedEle = [];
uDIonName = [];
uSpeciesM = [];
cnt = 1;
for i = 1:size(D, 2)
    idx = find( D(:, i) ~=0 );
    uDTemplate( idx, i ) = 1;
    oldIdx = find( DTemplate(:, i ) == 1 );
    ins = ismember( oldIdx, idx );
    if ~isempty(find(ins==1, 1))
        uDIonName{cnt} = DIonName{i}(ins, :);
        uSpeciesM(cnt) = SpeciesM(i);
        cnt = cnt + 1;
    end
    if ~isempty(find(ins==1, 1))
        usedEle = [usedEle i];
    end
end
uDTemplate = uDTemplate(:, usedEle);
uD = D(:, usedEle);

W = expRec.W(usedEle, :);
W(W<1e-2) = 0;
uWInfo = cell( size(W(:, :), 2 ), 1 );
dfW = zeros( length(uWInfo), 1 );
for i = 1:length(uWInfo)
    uWInfo{i} = find(W(:, i) ~=0);
    dfW(i) = length(uWInfo{i});
end


end

