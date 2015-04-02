function [ uDTemplate, uDIonName, uSpeciesM ] = improveDTemplate( mode, DTemplate, DIonName, speciesM )
%--------------------------------------------------------------------------
% modifyDTemplate: modify DTemplate
%--------------------------------------------------------------------------
% DESCRIPTION:
%   Modify dictionary template based on some chemcial knowledge. For
%   negative mode, the rule is [M+Na-2H], [M+H2O-H],
%   and [M+K-2H] must always appear with [M-H]
%   For positive mode, the rule is [M+2Na-H] must appear with [M+H] and
%   [M+Na], and [M+2K-H] must appear with [M+H] and [M+K]
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
        for j = 1:size(cName, 2)
            cc = cName(:,j);
            hFlag = detectPatInCell( cc, 'M\+H' );
            hNaFlag = detectPatInCell( cc, 'M\+Na' );
            if (hFlag && hNaFlag) == 0
                hNa2Flag = detectPatInCell( cc, 'M\+2Na-H' );
                if hNa2Flag == 1
                    insVec(j) = 1;
                end
            end
            hKFlag = detectPatInCell( cc, 'M\+K' );
            if (hFlag && hKFlag) == 0
                hK2Flag = detectPatInCell( cc, 'M\+2K-H' );
                if hK2Flag == 1
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
end

uDTemplate = DTemplate(:, keptVec);
uDIonName = DIonName(keptVec);
uSpeciesM = speciesM(keptVec);
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
