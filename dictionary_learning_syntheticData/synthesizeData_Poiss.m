function [simData] = synthesizeData_Poiss( SLEN, MLEN, HEIGHT, WIDTH, DTemplate, options, verbose )
%synthesizeData_Poiss synthesize data with poiss noise
%SLEN, MLEN, HEIGHT, WIDTH, the dimension of the data
%DTemplate, the dictionary template from which we generate gD (a zero-one
%matrix)
%verbose, 1 means show processing messages
%options: a data structure with the folowing possible fields:
%   WType: type of w, that is one of 'random', 'diffusion', 'thresholding',
%   default is 'diffusion'
%   WScale: scale of w, default is 1
%   supNum: the maximun # support, the maximun
%   # dictionary elements used in each sample, default is MLEN
%   coheMax: the maximun value of coherence of each dictionary elements,
%   default is 1
%   nonZMax: the percentage of nonZ entry in each dictionary elements, default
%   is 1
%verbose, 1 means show processing messages

%% setting input parameters
if ~isfield( options, 'WType' )
    WType = 'diffusion';
else
    WType = options.WType;
end
if ~isfield( options, 'WScale' )
    WScale = 1;
else
    WScale = options.WScale;
end
if ~isfield( options, 'supNum' )
    supNum = MLEN;
else
    supNum = options.supNum;
end
if ~isfield( options, 'coheMax' )
    coheMax = 1;
else
    coheMax = options.coheMax;
end
if ~isfield( options, 'nonZMax' )
    nonZMax = 1;
else
    nonZMax = options.nonZMax;
end

%% Generate orthogonal dictionary (random normal dist. generation)
fprintf( 'Synthesize D...\n' );
if ~isempty(DTemplate)
    [ gD, uDTemplate, usedTerm ] = synthesizeD( DTemplate, nonZMax, coheMax, supNum, verbose );
else
    [ gD, uDTemplate ] = synthesizeD_woTem( SLEN, MLEN, supNum, coheMax, verbose );
end

%% Generate w
fprintf( 'Synthesize W...\n' );
[ gW, gW0 ] = synthesizeW( SLEN, MLEN, HEIGHT, WIDTH, WType, WScale, verbose);

%% Generate Y
gY = genY( gD, gW, gW0 );
% gY = round(lambda);
simData.gY = gY; simData.gW = gW; simData.gW0 = gW0; simData.gD = gD; simData.uDTemplate = uDTemplate;
if DTemFlag == 1
    simData.usedTerm =usedTerm;
end
end
