function [simData] = synthesizeData_Poisson( SLEN, MLEN, HEIGHT, WIDTH, DTemplate, DOptions, WOptions, verbose )
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
params.sparsePrec = sparsePrec;
params.coheMax = coheMax;
if ~isempty(DTemplate)
    [SLEN, MLEN] = size(DTemplate);
    params.DTemplate = DTemplate;
    [ gD, uDTemplate ] = synthesizeD( 'DTemplate', params, verbose );
else
    params.SLEN = SLEN; params.MLEN = MLEN;
    [ gD, uDTemplate ] = synthesizeD( 'random', params, verbose );
end

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
if ~isfield( WOptions, 'supPerc' )
    supPerc = 0.2;
else
    supPerc = WOptions.supPerc;
end
if ~isfield( WOptions, 'sparsePerc' )
    if strcmp( type, 'random' ) == 0
        sparsePerc = 1;
    else
        sparsePerc = 0.2;
    end
else
    sparsePerc = WOptions.sparsePerc;
end

%% Generate w
fprintf( 'Synthesize W...\n' );
[ gW, gW0, usedElement ] = synthesizeW( MLEN, HEIGHT, WIDTH, type, supPerc, sparsePerc, scale, verbose);


%% Generate Y
gY = genY_Poisson( gD, gW, gW0 );
% gY = round(lambda);
simData.gY = gY; simData.gW = gW; simData.gW0 = gW0; simData.gD = gD; simData.uDTemplate = uDTemplate; simData.usedElement = usedElement;
end
