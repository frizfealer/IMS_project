function [simData] = synthesizeData_Poisson( LINK_FUNC, CONSTRAINTS, SLEN, MLEN, HEIGHT, WIDTH, DTemplate, DOptions, WOptions, verbose )
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

%% setting D's input parameters
if ~isfield( DOptions, 'sparsePrec' )
    if ~isempty(DTemplate)
        sparsePrec = 0.5;
    else
        sparsePrec = 0.1;
    end
else
    sparsePrec = DOptions.sparsePrec;
end
if ~isfield( DOptions, 'coheMax' )
    coheMax = 0.6;
else
    coheMax = DOptions.coheMax;
end

%% Generate orthogonal dictionary (random normal dist. generation)
fprintf( 'Synthesize D...\n' );
params.sparsePrec = sparsePrec;
params.coheMax = coheMax;
if ~isempty(DTemplate)
    [SLEN, MLEN] = size(DTemplate);
    params.DTemplate = DTemplate;
    [ gD, uDTemplate ] = synthesizeD( CONSTRAINTS, 'DTemplate', params, verbose );
else
    params.SLEN = SLEN; params.MLEN = MLEN;
    [ gD, uDTemplate ] = synthesizeD( CONSTRAINTS, 'random', params, verbose );
end
MLEN = size(gD, 2);
%% setting W'sinput parameters
if ~isfield( WOptions, 'type' )
    type = 'diffusion';
else
    type = WOptions.type;
end
if ~isfield( WOptions, 'scale' )
    scale = 1;
else
    scale = WOptions.scale;
end
if ~isfield( WOptions, 'supPrec' )
    supPrec = 0.2;
else
    supPrec = WOptions.supPrec;
end
if ~isfield( WOptions, 'sparsePrec' )
    if strcmp( type, 'random' ) == 0
        sparsePrec = 1;
    else
        sparsePrec = 0.2;
    end
else
    sparsePrec = WOptions.sparsePrec;
end

%% Generate w
fprintf( 'Synthesize W...\n' );
[ gW, gW0, usedElement ] = synthesizeW( MLEN, HEIGHT, WIDTH, type, supPrec, sparsePrec, scale, verbose);


%% Generate Y
gY = genY_differentNoise( LINK_FUNC, gD, gW, gW0 );
% gY = round(lambda);
simData.gY = gY; simData.gW = gW; simData.gW0 = gW0; simData.gD = gD; simData.uDTemplate = uDTemplate; simData.usedElement = usedElement;
end
