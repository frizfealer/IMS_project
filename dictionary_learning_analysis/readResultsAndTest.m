function [ lVec, tVec, pVec, valVec ] = readResultsAndTest( resFolder, inputPath )
%readResultsAndTest Summary of this function goes here
%   Detailed explanation goes here
load( inputPath );
cd( resFolder );
listing = dir( resFolder );
lVec = zeros( length(listing) - 2, 1 );
tVec = zeros( length(listing) - 2, 1 );
pVec = zeros( length(listing) - 2, 1 );
valVec = zeros( length(listing) - 2, 1 );
for i = 3:length(listing)
    fName = listing(i).name;
    load(fName);
    lVec(i-2) = expRec.lambda;
    tVec(i-2) = expRec.theta;
    pVec(i-2) = expRec.phi;
    [ valVec(i-2) ] = validationOnTesting( aMatrix, dataCube, expRec.D );
end

end

