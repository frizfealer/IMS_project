function [simData] = synthesizeData_DPatterns_Poiss( SLEN, MLEN, HEIGHT, WIDTH, TemVar, nTemVar, type, scale, verbose )
%synthesizeData_Poiss synthesize data with poiss noise
%SLEN, MLEN, HEIGHT, WIDTH, the dimension of the data
%DTemplate, the dictionary template from which we generate gD
%sparseVal, the sparse value of dictionary entries, 0~1, 1 means no sparse,
%at least used one signal in an element, default is 1
%coheVal, the maximun allowed coherence values of dictionary elements,
%default is 1
%type, type of w, 'random', 'diffusion', 'thresholding'
%scale, the scale of w, default is 1
%verbose, 1 means show processing messages
if isempty( scale )
    scale = 1;
end

DTemFlag = 0;
if isempty(TemVar)
    supNum = nTemVar.supNum;
    if isempty( nTemVar.coheVal )
        coheVal = 1;
    else
        coheVal = nTemVar.coheVal;
    end
else
    DTemplate = TemVar.DTemplate;
    sparseVal = TemVar.sparseVal;
    if isempty( TemVar.coheVal )
        coheVal = 1;
    else
        coheVal = TemVar.coheVal;
    end
    DTemFlag = 1;
end
if DTemFlag == 1
    assert( SLEN == size(DTemplate, 1) );
    assert( MLEN <= size(DTemplate, 2) );
end
%% Generate w
fprintf( 'Synthesize W...\n' );
[ gW, gW0 ] = synthesizeW( SLEN, MLEN, HEIGHT, WIDTH, type, scale, verbose);
%% Generate orthogonal dictionary (random normal dist. generation)
fprintf( 'Synthesize D...\n' );
if DTemFlag == 1
    [ gD, uDTemplate, usedTerm ] = synthesizeD( SLEN, MLEN, DTemplate, sparseVal, coheVal, verbose );
else
    [ gD, uDTemplate ] = synthesizeD_woTem( SLEN, MLEN, supNum, coheVal, verbose );
end
%% Generate Y
gY = genY( gD, gW, gW0 );
% gY = round(lambda);
simData.gY = gY; simData.gW = gW; simData.gW0 = gW0; simData.gD = gD; simData.uDTemplate = uDTemplate;
if DTemFlag == 1
    simData.usedTerm =usedTerm;
end
end
