function [ rSet ] = estimateDRelevant( dataCube, outD, DTemplate, indMap, percent, intThres )
%estimateRelevant estimate relevant of each dictionary elements
%dataCube: input MASS signals
%outD: current dictionary
%DTemplate: dictionary template (0-1 entry)
%indMap: indice map of the same size of sample, is a data field of BlkDS
%percent: the percent of dictionary elements you want to use in all outer
%iterations
% intThres: the percentage of intensity you want to include in your m/z
% selections
if isempty( percent )
    percent = 0.2;
end
if isempty( intThres )
    intThres = 0.8;
end
mLen = size( outD, 2 );
[sLen, nLen] = size(dataCube);
% opts=[];
% opts.init=2; % starting from a zero point
% opts.tFlag=5;       
% opts.maxIter=100; 
% opts.nFlag=0;       % without normalization
% opts.rFlag=1;       % the input parameter 'rho' is a ratio in (0, 1)
% %opts.rsL2=0.01;     % the squared two norm term
% 
% %----------------------- Run the code LeastR -----------------------
% fprintf('\n mFlag=0, lFlag=0 \n');
% opts.mFlag=0;       % treating it as compositive function 
% opts.lFlag=0;       % Nemirovski's line search

fprintf( 'computing residauls of each dictionary elements...\n' );
Y = dataCube(:, :);
Y = log(Y); Y(Y==-inf)=0;
resAry = zeros( mLen, 1 );
W = zeros( nLen, 1 );
for i = 1:mLen
    cIdx = find( DTemplate(:,i) ==1 );
    cY = Y(cIdx, :);
    cDEle = outD(cIdx, i);
    parfor j = 1:nLen
        if indMap(j) == 1
            if ~isempty( find( cY(:,j) ~= 0, 1 ) ) 
%         [W(1,j), ~, ~]= LeastR(outD(cIdx,i), Y(cIdx, j), 0, opts);
                W(j) = cY(:, j) \ cDEle;
            end
        end
    end
%     if W(1,j) < 0
%         keyboard();
%     end
    res = cY -cDEle*W';
    res = res(:);
    resAry(i) = norm(res);
end

fprintf( 'constructing possible dictionary elements of each m/z...\n' );
ins = zeros( sLen, 1 );
for i = 1:sLen
    ins(i) = sum( dataCube(i, :) );
end
[sIns,idx] = sort( ins, 'descend' );
totalIonInt = sum(ins);
cSum = 0;
for i = 1:sLen
    cSum = cSum + sIns(i);
    if cSum / totalIonInt > intThres
        break;
    end
end
if i == 1
    intEleVec = idx(1);
else
    intEleVec = idx(1:(i-1));
end
pEleList = [];
for i = 1:length(intEleVec)
    pEle = find( DTemplate( intEleVec(i), : ) == 1 );
    cEleRes = resAry(pEle);
    [~,idx] = sort(cEleRes);
    pEleList{i} = pEle(idx);
end

%round-robin manner to choose dictionary elements
fprintf( 'round-robinly choose possible dictionary elements...\n' );
eleLimit = round( mLen*percent );
rSet = [];
rRecVec = ones( length(intEleVec), 1 );
while length(rSet) < eleLimit
    for i = 1:length(intEleVec)
        cPEle = pEleList{i};
        if rRecVec(i) <= length(cPEle)
            rSet = unique( [rSet; cPEle(rRecVec(i))] );
            rRecVec(i) = rRecVec(i) + 1;
        end
    end
end

end

