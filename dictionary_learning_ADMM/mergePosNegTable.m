function [ tSpeciesM, tDTemplate, tDIonName ] = mergePosNegTable( pSpeciesM, nSpeciesM, pDTemplate, nDTemplate, pDIonName, nDIonName, errRange, MPPATH, inpeakMZ )
%mergePosNegTable merge Positive and Negative DTemplate tables
%parameters:
%pSpeciesM, nSpeciesM, pDTemplate: output from positive table.
%nDTemplate, pDIonName, nDIonName: output from negative table.
%errRange. the tolerance of m/z that to consider species MW is the same
%MPPATH: the file path of negative ion mode.
%inpeakMZ: the m/z axis of negative mode data.

fileID = fopen(MPPATH);
C = textscan(fileID, '%f %s', 'delimiter', {','});
fclose(fileID);
mpAry = C{1,1};
nameAry = C{1,2};

tSpeciesM = pSpeciesM;
tDTemplate = zeros( size( pDTemplate, 1 ) + size( nDTemplate, 1 ), size( pDTemplate, 2 ) );
tDTemplate(1:size(pDTemplate,1),:) = pDTemplate;
tDIonName = pDIonName;
pLen = length( tDIonName );
cnt = 1;
for i = 1:length(nSpeciesM)
    tmp = abs( pSpeciesM - nSpeciesM(i) );
    [res, idx] = min(tmp);
    if res < errRange
        possibleRelatedPeaks = mpAry+pSpeciesM(idx);
        [ mVec ] = computeMapping( possibleRelatedPeaks, inpeakMZ, errRange );
        curPeak = mVec( ~isnan( mVec(:, 2) ), 2 );
        tDTemplate(curPeak+size(pDTemplate, 1), idx) = 1;
        cName = tDIonName{idx};
        tDIonName{idx} = [cName; nameAry( ~isnan( mVec(:,2) ) )];
    else
        tSpeciesM = [tSpeciesM nSpeciesM(i)];
        tDTemplate = [tDTemplate [zeros( size(pDTemplate, 1), 1); nDTemplate(:,i)] ];
        tDIonName{pLen+cnt} = nDIonName{i};
        cnt = cnt + 1;
    end
end

fprintf( 'Computing redundancy in the DTemplate...\n' );
%redundent is defined as:
% A:[1 1 0 0]
% B:[1 1 1 0]
% A is redudent of B.
redunVec = zeros( size( tDTemplate, 2 ), 1 );
% recVec = zeros( size(DTemplate, 2), 1);
for i = 1:size(tDTemplate, 2)
    cSet = tDTemplate(:,i);
    for j = i:size(tDTemplate, 2)
        tmp = cSet - tDTemplate(:,j);
        if max(tmp)-min(tmp)==1 && min(tmp) == -1
            redunVec(i) = 1;
%             recVec(i) = j;
            break;
        end
    end
end
tmp = tDIonName;
tDIonName = {};
cnt = 1;
for i = 1:length(tmp)
    if redunVec(i) == 0
        tDIonName{cnt} = tmp{i};
        cnt = cnt + 1;
    end
end
tDTemplate = tDTemplate(:, redunVec==0);
tSpeciesM = tSpeciesM(redunVec == 0);

end

