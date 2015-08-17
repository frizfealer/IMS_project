function [ DL ] = initDL( WA, DI )
%--------------------------------------------------------------------------
%initDInterior: initialize the data structure of Dictionary Learning
%--------------------------------------------------------------------------
% DESCRIPTION:
%
% INPUT ARGUMENTS:
%   WA, the data structue of W ADMM
%   DI, the data structure of D Interior
% OUTPUT ARGUMENTS:
%   DL, the data structure of dictionary learning, has the following
%   fields,
%       WA, the W ADMM.
%       DI, the D Interiror.
%       itNum: # loop of the dictionary learning (outer loop).
%       SAVE_TEM_PERIOD, how many iterations to save current model.
%       LP_TOL, a termination criterion. If the difference of LP between
%       the recent two iterations is <= current LP * LP_TOL, then break.
%       W_TOL, D_TOL, other termination criteria. If the difference of D is
%       <= D_TOL and the difference of W is <= max(current W)*W_TOL, then
%       break.
%       D_HIST_FLAG, a flag indicates whether to save D in each iteration.
%       W_HIST_FLAG, a flag indicates whether to save W in each iteration.
%       LPAry, a vector recording log posterior values for each iteration.
DL.W = WA;
DL.D = DI;
DL.itNum= 100;
DL.SAVE_TEM_PERIOD = 3;
DL.LP_TOL = 1e-6;
DL.W_TOL = 1e-4;
DL.D_TOL = 1e-4;
DL.D_HIST_FLAG = 0;
DL.W_HIST_FLAG = 0;
DL.LPAry = zeros( 1, DL.itNum+1 );
end

