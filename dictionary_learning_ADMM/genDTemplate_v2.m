function [ DTemplate, DIonName, speciesM ] = genDTemplate_v2( dataCube, inpeakMZ, MPPATH, smallestM, errRange, varargin )
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

interval = mpAry - mpAry(1);
mapping = zeros( length(interval), length(interval) );
for i = 2:length(interval)
    mapping(i,:) = interval - interval(i);
end
mapping(1,:) = interval;

ins = zeros( size(dataCube, 1), 1 );
for i = 1:size(ins)
    ins(i) = sum( dataCube(i, :) );
end
[~, ord] = sort(ins, 'descend');
oInMolMZ = inMolMZ(ord);
speciesM = [];
DTemplate = [];
DIonName = {};
cnt = 1;
fprintf( 'generating DTemplate...\n' );
while ~isempty(oInMolMZ)
    cMZ = oInMolMZ(1);
    for j = 1:size(mapping, 1)
         possibleRelatedPeaks = mapping(j,:)+cMZ;
            [ mVec ] = computeMapping( possibleRelatedPeaks', inpeakMZ', errRange );
            curPeak = mVec( ~isnan( mVec(:, 2) ), 2 );
            if isempty( DTemplate )
                DTemplate = zeros( SLEN, 1 );
                DTemplate(curPeak, 1) = 1;
                DIonName{cnt} = nameAry( ~isnan( mVec(:,2) ) );
                speciesM = [speciesM cMZ-mpAry(j)];
                
                idx =  ismember(oInMolMZ, inMolMZ(curPeak)) ;
                oInMolMZ(idx) = [];
            else
                curTem = zeros( SLEN, 1 );
                curTem( curPeak, 1 ) = 1;
                inFlag = 0;
                for k = size( DTemplate, 2 ):-1:1
                    tmp = DTemplate(:, k) - curTem;
                    ins = sum(abs(tmp));
                    if ins ==0
                        inFlag = 1;
                        break;
                    end
                end
                if inFlag == 0
                    cnt = cnt + 1;
                    DTemplate = [DTemplate curTem];
                    DIonName{cnt} = nameAry( ~isnan( mVec(:,2) ) );
                    speciesM = [speciesM cMZ-mpAry(j)];
                    idx =  ismember(oInMolMZ, inMolMZ(curPeak)) ;
                    oInMolMZ(idx) = [];
                end
            end
    end
end

end
%
% testing
% for i = 1:length(DIonName)
%     if size(DIonName{i}, 1) ~= length(find(DTemplate(:,i)==1))
%         fprintf( '%d\n', i );
%     end
% end
