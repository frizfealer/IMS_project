function [ uDTemplate, uDIonName, uSpeciesM ] = improveDTemplate( mode, DTemplate, DIonName, speciesM )
%--------------------------------------------------------------------------
% modifyDTemplate: modify DTemplate
%--------------------------------------------------------------------------
% DESCRIPTION:
%   Modify dictionary template based on some chemcial knowledge. For
%   negative mode, the rule is if a dictionary element has one of the three 
%   ions( [M+Na-2H], [M+H2O-H], and [M+K-2H]) but does not have [M-H], then 
%   we remove this dictionary element.
%   For positive mode, the rule is 
%   1. if a dictionary element has [M+2Na-H] and does not have [M+H] and 
%   [M+Na], delete this element.
%   2. if a dictionary element has [M+2Na-H] and have one of the two ions  
%   ([M+H] and [M+Na]), retain this element but delete [M+2Na-H].
%   3. if a dictionary element has [M+2K-H] and does not have [M+H] and 
%   [M+K], delete this element.
%   4. if a dictionary element has [M+2K-H] and have one of the two ions  
%   ([M+H] and [M+K]), retain this element but delete [M+2K-H].
%
% INPUT ARGUMENTS:
%   mode, positve mode or negative mode
%   DTemplate, DIonName, speciesM, the dictionary template variables
%
% OUTPUT ARGUMENTS:
%   uDTemplate, uDIonName, uSpeciesM, the updated variables

keptVec = [];
if strcmp( mode, 'negative') == 1
    for i = 1:length(DIonName)
        cName = DIonName{i};
        % a dictioanry element may have more than one ion form
        % combinations, need to test all of them
        insVec = zeros( size(cName, 2), 1 );
        for j = 1:size(cName, 2)
            cc = cName(:,j);
            hFlag = detectPatInCell( cc, 'M-H' );
            if hFlag == 0
                hNaFlag = detectPatInCell( cc, 'M\+Na-2H' );
                hKFlag = detectPatInCell( cc, 'M\+K-2H' );
                hHFlag = detectPatInCell( cc, 'M-H2O-H' );
                if hNaFlag == 1 || hKFlag == 1 || hHFlag == 1
                    insVec(j) = 1;
                end
            end
        end
        oFlag = insVec(1);
        for j = 2:size(cName, 2)
            oFlag = oFlag *insVec(j);
        end
        if oFlag == 0
            keptVec = [keptVec i];
        end
        DIonName{i}(:, insVec==1) = [];
        speciesM{i}(:, insVec==1) = [];
    end
elseif strcmp( mode, 'positive' ) == 1
    for i = 1:length(DIonName)
        cName = DIonName{i};
        insVec = zeros( size(cName, 2), 1 );
        % a dictioanry element may have more than one ion form
        % combinations, need to test all of them
        for j = 1:size(cName, 2)
            cc = cName(:,j);
            hFlag = detectPatInCell( cc, 'M\+H' );
            hNaFlag = detectPatInCell( cc, 'M\+Na' );
            if (hFlag == 0 && hNaFlag == 0)
                hNa2Flag = detectPatInCell( cc, 'M\+2Na-H' );
                if hNa2Flag == 1
                    insVec(j) = 1;
                end
            elseif (hFlag == 1 && hNaFlag == 0) || (hFlag == 0 && hNaFlag == 1)
                hNa2Flag = detectPatInCell( cc, 'M\+2Na-H' );
                if hNa2Flag == 1
                    insVec(j) = 2;
                end
            end
            hKFlag = detectPatInCell( cc, 'M\+K' );
            if (hFlag == 0 && hKFlag == 0)
                hK2Flag = detectPatInCell( cc, 'M\+2K-H' );
                if hK2Flag == 1
                    insVec(j) = 1;
                end
            elseif (hFlag == 1 && hKFlag == 0) || (hFlag == 0 && hKFlag == 1)
                hK2Flag = detectPatInCell( cc, 'M\+2K-H' );
                if hK2Flag == 1
                    insVec(j) = 2;
                end
            end
        end
        oFlag = insVec(1);
        for j = 2:size(cName, 2)
            oFlag = oFlag *insVec(j);
        end
        if oFlag == 0 %one of the combination is correct
            keptVec = [keptVec i];
            %remove all other possible wrong combination
            DIonName{i}(:, insVec==2) = [];
            speciesM{i}(:, insVec==2) = [];
            DIonName{i}(:, insVec==1) = [];
            speciesM{i}(:, insVec==1) = [];
        elseif mod(oFlag, 2) == 0 
            %one of the combination is partially corret
            %keep this element.
            keptVec = [keptVec, i];
            %find the first correct combination.
            %Actually in current setting, every combinations in this case
            %is the same.
            conIdx = find( insVec == 2, 1);
            %get the ion forms in this combination.
            cName = DIonName{i}(:, conIdx);
            %detect M+2Na-H position.
            NaPos = -1;
            for j = 1:length(cName)
                cc = cName(j);
                f1 = detectPatInCell( cc, 'M\+2Na-H' );
                if f1 == 1
                    NaPos = j;
                    break;
                end
            end
            %detect M+2K-H position.
            KPos = -1;
            for j = 1:length(cName)
                cc = cName(j);
                f1 = detectPatInCell( cc, 'M\+2K-H' );
                if f1 == 1
                    KPos = j;
                    break;
                end
            end
            %remove incorrect ion forms in the dictionary template
            idx = find(DTemplate(:, i)==1);
            cSet = 1:length(DIonName{i}(:, conIdx));
            if NaPos ~= -1
                DTemplate( idx(NaPos), i ) = 0;
            end
            if KPos ~= -1
                DTemplate( idx(KPos), i ) = 0;
            end
            %only retain one correct ion forms 
            %in the combination name lists.
            cSet = setdiff( cSet, NaPos );
            cSet = setdiff( cSet, KPos );
            uCon = DIonName{i}(cSet, conIdx);
            DIonName{i}= uCon;
            %keep all possible putative molecular weight
            %not sure if it is reasonable.
            speciesM{i} = speciesM{i}(conIdx);
        end
    end
end

uDTemplate = DTemplate(:, keptVec);
uDIonName = DIonName(keptVec);
uSpeciesM = speciesM(keptVec);

%remove redunduct dictionary elements
rmVec = [];
for i = 1:size(uDTemplate, 2)-1
    d1 = uDTemplate(:, i);
    for j = i+1:size(uDTemplate, 2)
        d2 = uDTemplate(:, j);
        tmp = d1 - d2;
        ins = sum(abs(tmp));
        if ins == 0
            rmVec = [rmVec, j];
            uDIonName{i} = [uDIonName{i} uDIonName{j}];
            uSpeciesM{i} = [uSpeciesM{i} uSpeciesM{j}];
        end
    end
end
uDTemplate(:, rmVec) = [];
uDIonName(rmVec) = [];
uSpeciesM(rmVec) = [];

mvVec = zeros( length(uSpeciesM), 1 );
for i = 1:length(uSpeciesM)
    mvVec(i) = mean(uSpeciesM{i});
end
[~, idx] = sort( mvVec );
uSpeciesM = uSpeciesM(idx);
uDIonName = uDIonName(idx);
uDTemplate = uDTemplate(:, idx);

end

function [flag] = detectPatInCell( cc, pat )
flag = 0;
idx = regexp( cc, ['^', pat, '$'] );
for z = 1:length(idx)
    if ~isempty(idx{z})
        flag = 1;
        break;
    end
end
end
