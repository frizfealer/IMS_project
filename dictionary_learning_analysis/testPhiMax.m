function [ dfVec, phiVec ] = testPhiMax( Y, W, W0, D, DTemplate, startVal, intVal, intThres, testNum, seq )
%testLambdaMax test lambda max
if ~isempty(seq)
    phiVec = seq;
else
    endVal = startVal+testNum*intVal;
    phiVec = startVal:intVal:endVal;
end
dfVec = zeros(length(phiVec), 1);
parfor i = 1:length(phiVec)
    phi = phiVec(i);
    fprintf( 'phi = %g\n', phi );
    [~, hei, wid] = size(Y);
    aMatrix = ones( hei, wid );
    itNum = 200;
    scaleFactor = computeScaleFactor( Y, aMatrix );
    BlkDS = conBLKDS( Y );
    validMap = BlkDS.indMap .* aMatrix;
    [ uD ]= updateD_v8_ipopt( Y, W, W0, D, DTemplate, validMap, 1, phi, scaleFactor, itNum, 1e-2 );    
    dfVec(i) = length(find(uD(:)>intThres));
end

end

