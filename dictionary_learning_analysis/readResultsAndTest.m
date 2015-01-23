function [ lVec, tVec, pVec, valVec ] = readResultsAndTest( resFolder, inputPath, linkFunc, parForFlag )
%readResultsAndTest Summary of this function goes here
%   Detailed explanation goes here
load( inputPath );
cd( resFolder );
listing = dir( pwd );
lVec = zeros( length(listing) - 2, 1 );
tVec = zeros( length(listing) - 2, 1 );
pVec = zeros( length(listing) - 2, 1 );
valVec = zeros( length(listing) - 2, 1 );
expRecCell = cell( length(listing)-2, 1 );
fprintf( 'loading all the results...' );
for i = 3:length(listing)
    fName = listing(i).name;
    load(fName);
    expRecCell{i-2} = expRec;
end
if parForFlag == 1
    parfor i = 1:length(expRecCell)
        fprintf('%d\n', i);
        expRec = expRecCell{i};
        lVec(i) = expRec.lambda;
        tVec(i) = expRec.theta;
        pVec(i) = expRec.phi;
        [ valVec(i) ] = validationOnTesting( aMatrix, dataCube, expRec.D, expRec.W, expRec.W0, lVec(i), linkFunc );
    end
elseif parForFlag == 0
    for i = 1:length(expRecCell)
        fprintf('%d\n', i);
        expRec = expRecCell{i};
        lVec(i) = expRec.lambda;
        tVec(i) = expRec.theta;
        pVec(i) = expRec.phi;
        [ valVec(i) ] = validationOnTesting( aMatrix, dataCube, expRec.D, expRec.W, expRec.W0, lVec(i), linkFunc );
    end
end

end

