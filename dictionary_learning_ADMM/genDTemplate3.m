function [ DTemplate, DIonName, speciesM ] = genDTemplate3( inpeakMZ, MPPATH, smallestM, errRange, usedIonForms, varargin )
%genDTemplate generate DTemplate within epsilon error
%inpeakMZ [sLen 1],
%MPPATH: the pattern file path,
%'D:\RA\Project_dictionaryLearning_IMS\molecule_profile_new_trim.csv'
%smallestM: the smallest molecule's weight
%errRange: the error tolerance to bin m/z to inpeakMZ
%inMolMZ: the additional parameter that users can specified the m/z to
%generate peaks, instead of using the whole inpeakMZ
%return DTemplate, a template with 1 means possible and 0 means impossible
%DIonName, "possible" adducts generated that signals. This information is
%just one of the possible answers.
fileID = fopen(MPPATH);
C = textscan(fileID, '%f %s', 'delimiter', {','});
fclose(fileID);
mpAry = C{1,1};
nameAry = C{1,2};
SLEN = length(inpeakMZ);
if length(varargin) == 1
    inMolMZ = varargin{1};
else
    inMolMZ = inpeakMZ;
end
assert( max(usedIonForms) < length(mpAry) );




interval = mpAry - mpAry(1);
mapping = zeros( length(usedIonForms), length(interval) );
for i = 1:length(usedIonForms)
    mapping(i,:) = interval - interval(usedIonForms(i));
end

speciesM = [];
DTemplate = [];
DIonName = {};
cnt = 1;
fprintf( 'generating DTemplate...\n' );
for i = 1:length(inpeakMZ)
%     fprintf( 'm/z %g\n', mMVec(i));
    curMZ = inpeakMZ(i);
    for j = 1:length(usedIonForms)
        cMapping = mapping(j, :);
         possibleRelatedPeaks = cMapping+curMZ;
         [ mVec ] = computeMapping( possibleRelatedPeaks', inpeakMZ, errRange );
         curPeak = mVec( ~isnan( mVec(:, 2) ), 2 );
         if isempty( DTemplate )
             DTemplate = zeros( SLEN, 1 );
             DTemplate(curPeak, 1) = 1;
             DIonName{cnt} = nameAry( ~isnan( mVec(:,2) ) );
             speciesM{cnt} = (possibleRelatedPeaks(1)-mpAry(1));
         else
             curTem = zeros( SLEN, 1 );
             curTem( curPeak, 1 ) = 1;
             inFlag = 0;
             for k = size( DTemplate, 2 ):-1:1
                 tmp = DTemplate(:, k) - curTem;
                 ins = sum(abs(tmp));
                 if ins ==0
                     inFlag = 1;
                     speciesM{k} = [speciesM{k} possibleRelatedPeaks(1)-mpAry(1)];
                     DIonName{k} = [DIonName{k}  nameAry( ~isnan( mVec(:,2) ) )];
                 end
             end
             if inFlag == 0
                 cnt = cnt + 1;
                 DTemplate = [DTemplate curTem];
                 DIonName{cnt} = nameAry( ~isnan( mVec(:,2) ) );
                 speciesM{cnt} = (possibleRelatedPeaks(1)-mpAry(1));
             end
         end
    end
end

fprintf( 'computing redundancy in the DTemplate...\n' );
%redundent is defined as:
% A:[1 1 0 0]
% B:[1 1 1 0]
% A is redudent of B.
redunVec = zeros( size( DTemplate, 2 ), 1 );
% recVec = zeros( size(DTemplate, 2), 1);
nonZeroEnt = cell( size(DTemplate, 2), 1 );
nonZeroLen = zeros( size(DTemplate, 2), 1 );
for i = 1:size(DTemplate, 2)
    nonZeroEnt{i} = find( DTemplate(:, i) == 1 );
    nonZeroLen(i) = length(nonZeroEnt{i});
end

for i = 1:size(DTemplate, 2)
    fprintf('%d\n', i );
    cNZLen = nonZeroLen(i);
    cEle = DTemplate(:, i);
    targetEle = find( nonZeroLen > cNZLen );
    for j = 1:length(targetEle)
        oEle = DTemplate(:, targetEle(j));
        ins = cEle | oEle;
        if length(find(ins==1)) == nonZeroLen(targetEle(j))
           redunVec(i) = 1;
            break;
        end
    end
end
tmp1 = DIonName;
tmp2 = speciesM;
DIonName = {};
speciesM = {};
cnt = 1;
for i = 1:length(tmp1)
    if redunVec(i) == 0
        DIonName{cnt} = tmp1{i};
        speciesM{cnt} = tmp2{i};
        cnt = cnt + 1;
    end
end
DTemplate = DTemplate(:, redunVec==0);
end

