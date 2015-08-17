function [ DInterior_class ] = initDInterior( DTemplate, dataCube, DIonName, SpeciesM, init_method, CONSTRAINT )
%--------------------------------------------------------------------------
%initDInterior: initialize the data structure of D Interior point
%--------------------------------------------------------------------------
% DESCRIPTION:
%
% INPUT ARGUMENTS:
%   DTemplate, dictionary template.
%   DIonName, the ion names for each element in the template.
%   SpeciesM, the putative MWs for each element in the template.
%   init_method, the method for dictionary initialization.
% OUTPUT ARGUMENTS:
%   DInterior_class, param data structre, has fields of:
%       itNum: # loop of ADMM update on W.
%       phi: a hyperparameter for controling D's sparsity.
%       DTemplate, DIonName, SpeciesM, from the inputs.
%       HESFlag, a flag indicate whether to use Hessian in optimization,
%       under construction.
%       W_LOWER_BOUND: the lower bound of the element' abundance to be
%       considered as a valid element. The max of the element should be
%       larger the max of the all element * W_LOWER_BOUND
%       LINK_FUNC: link function in the model. Default is 'identity'. 
%       Others are under construction.
%       kappa, parameter for the negative binomial distribution, under
%       construction.
%       CONSTRAINT, constraints on the dictionary elements, there are two
%       options: 'L2_SQUARE' and 'L1'
%       D, the dictionary
%       stuckFlag, a flag indicate whether the process stuck while
%       updating. "Stuck" means cannot make LP smaller within the
%       2*DInterior_class.itNum iterations.
DInterior_class.itNum= 200;
DInterior_class.phi = 0;
DInterior_class.DTemplate = DTemplate;
DInterior_class.DIonName = DIonName;
DInterior_class.SpeciesM = SpeciesM;
DInterior_class.HESFlag = 0;
DInterior_class.kappa = [];
DInterior_class.W_LOWER_BOUND = 1e-4;
DInterior_class.LINK_FUNC = 'identity';
if isempty(CONSTRAINT)
    DInterior_class.CONSTRAINT = 'L2_SQUARE';
else
    DInterior_class.CONSTRAINT = CONSTRAINT;
end
if isempty(init_method)
    DInterior_class.INIT_METHOD = 'NNMF';
else
    DInterior_class.INIT_METHOD = init_method;
end

DInterior_class.D = initD( dataCube, DTemplate, ...
    DInterior_class.INIT_METHOD, DIonName, DInterior_class.CONSTRAINT );

DInterior_class.stuckFlag = 0;

end

