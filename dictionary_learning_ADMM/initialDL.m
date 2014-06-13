function [ estD ] = initialDL( inY, tho, eps )
%initialDL Implementation of InitDictionaryLearn in paper
%"Learning Sparsely Used Overcomplete Dictionaries"
%tho shoud be 0.5 - s^2*mu0/d^0.5 > 0;
%eps shoud be 0.25 > eps^2 > 32sM^2/m^2 *...
iY = inY(:,:);
% corGraph = abs(iY'*iY);
corGraph = abs( corrcoef(iY) );
corGraph( corGraph >= tho ) = 1;
corGraph( corGraph ~= 1 ) = 0;
% estD = zeros(0, 0);
A = zeros(0, 0);
epsSq = (2*eps)^2;
cliqueGrp = {};
curNum = 1;
for i = 1:size(corGraph, 1)
    fprintf( '%d\n', i);
    for j = (i+1):size(corGraph, 1)
        if( corGraph( i, j ) == 1 )

            Ni = find( corGraph(i,:) == 1 );
            Nj = find( corGraph(j,:) == 1 );
            S = fastIntersect( Ni, Nj );
            if UniqueIntersect( S, corGraph )
                if isempty( cliqueGrp )
                    cliqueGrp{curNum} = S;
                    curNum = curNum + 1;
                else
                    inFlag = 0;
                    for k = length(cliqueGrp):-1:1
                        if length( cliqueGrp{k} ) == length(S) && norm(S-cliqueGrp{k},1) == 0
                            inFlag = 1;
                            break;
                        end
                    end
                    if inFlag == 0
                        cliqueGrp{curNum}=S;
                        curNum = curNum + 1;
                    end
                end
%                 Q = zeros( size(iY, 1), size(iY, 1) );
%                 for k = 1:length(S)
%                     Q = Q + iY(:,S(k))*iY(:,S(k))';
%                 end
%                 Q = iY(:,S)*iY(:,S)';
%                 tic
%                 u = pca(Q, 'NumComponents' , 1 );
%                 toc
%                 tic
                %[u,~,~] = svd(Q);
%                 [U,~,~]= svds(Q,1);
%                 [U,~] = nnmf(Q,1);
%                 toc
                
                %min|a-b|
%                 if isempty( A )
%                     A = U;
%                 else
%                     uMat = repmat(U, 1, size(A, 2) );
%                     
%                     if min( sum(( A - uMat ).^2, 1 ) ) > epsSq
%                         A = [A U];
%                     end
%                 end
            end
            
        end
    end
end
keyboard();
estD = A;

end

function [res] = UniqueIntersect(S, corG)
    thres = floor( length(S) / 2) * 61 / 64;
    res = 0;
    for i = 1:1e2
        edgeNum = 0;
%         ranVec = randperm(length(S));
        ranVec = RandInplacePermute( S );
        if mod( numel(ranVec), 2 ) ~= 0
            ranVec = ranVec(1:end-1);
        end
        ranVec = reshape( ranVec, numel(ranVec)/2, 2 );
        ranVec3 = ranVec(:,1) + ( ranVec(:,2) - 1 )*size(corG,1);
%         for j = 1:2:(length(ranVec) - 1)
%             edgeNum = edgeNum + corG( S(ranVec(j)), S(ranVec(j+1) ) );
%         end
%         for j = 1:length(ranVec1)-1
%             edgeNum = edgeNum + corG( S(ranVec2(j)), S(ranVec1(j) ) );
%         end
        edgeNum = sum( corG( ranVec3 ) );
        if edgeNum > thres
            res = 1;
            break;
        end
    end
end
