function [ lVec, tVec, pVec, valVec ] = readResultsAndTest( resFolder, inputPath, linkFunc, parForFlag )
%--------------------------------------------------------------------------
% readResultsAndTest: read results from all hyperparameters exp.
%--------------------------------------------------------------------------
% DESCRIPTION:
%   read the results from all hyperparameters experiments and tests
%	on the leave-out data.
%
% INPUT ARGUMENTS:
%   resFolder, result folder path
%	inputPath, basic variables path
%	linkFunc, should be 'identity', under construction
%	parForFlag for running regression
% OUTPUT ARGUMENTS:
%   lVec, tVec, pVec, valVec, the lambda, theta, phi, and their
%	log-posterior (the larger the better).
s = load(inputPath);
dataCube = s.IMSD.dataCube;
aMatrix = s.aMatrix;
listing = dir( resFolder );
lVec = zeros( length(listing) - 2, 1 );
tVec = zeros( length(listing) - 2, 1 );
pVec = zeros( length(listing) - 2, 1 );
valVec = zeros( length(listing) - 2, 1 );
expRecCell = cell( length(listing)-2, 1 );
fprintf( 'loading all the results...' );
for i = 3:length(listing)
    fName = listing(i).name;
    load([resFolder '/' fName]);
    expRecCell{i-2} = expRec;
end

if parForFlag == 1
    parfor i = 1:length(expRecCell)
        fprintf('%d\n', i);
        expRec = expRecCell{i};
        lVec(i) = expRec.lambda;
        tVec(i) = expRec.theta;
        pVec(i) = expRec.phi;
        if strcmp( linkFunc, 'negative_binomial' ) == 1
            [ valVec(i) ] = validationOnTesting( aMatrix, dataCube, expRec.D, expRec.W, expRec.W0, lVec(i), linkFunc, expRec.kappa );
        else
            [ valVec(i) ] = validationOnTesting( aMatrix, dataCube, expRec.D, expRec.W, expRec.W0, lVec(i), linkFunc );
        end
    end
elseif parForFlag == 0
    for i = 1:length(expRecCell)
        fprintf('%d\n', i);
        expRec = expRecCell{i};
        lVec(i) = expRec.lambda;
        tVec(i) = expRec.theta;
        pVec(i) = expRec.phi;
        if strcmp( linkFunc, 'negative_binomial' ) == 1
            [ valVec(i) ] = validationOnTesting( aMatrix, dataCube, expRec.D, expRec.W, expRec.W0, lVec(i), linkFunc, expRec.kappa );
        else
            [ valVec(i) ] = validationOnTesting( aMatrix, dataCube, expRec.D, expRec.W, expRec.W0, lVec(i), linkFunc );
        end
    end
end

end

