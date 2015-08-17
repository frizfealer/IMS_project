function [ InputData ] = initInputData( inData, mzAxis, leaveOutPercent )
%--------------------------------------------------------------------------
%initInputData: initialize the data structure of InputData
%--------------------------------------------------------------------------
% DESCRIPTION:
%
% INPUT ARGUMENTS:
%   inData, a three dimensional data cube
%   mzAxis, a vector that has m/z values correspond to inData.
%   leaveOutPercent, how many percentage in leave-out. It should be in
%   range of 0~100.
% OUTPUT ARGUMENTS:
%   InputData, a data structre, has fields of:
%       dataCube, the data cube.
%       BlkDS, the data structure of the position info. of this data.
%       Rblk, the fused relationship matrix for fused term W update.
%       scaleFactor, under construction, always set to 1.
%       logFY, log factorial of the data cube, used when computing log
%       posterior.
%       aMatrix, a 1-0 matrix that indicates whether to use a grid cell
%       when training. the cells that are one in the matrix is the subset
%       of the cells that are one in BlkDS.indMap

InputData.dataCube = inData;
InputData.mzAxis = mzAxis;
InputData.BlkDS = conBLKDS( inData );
[~, hei, wid] = size( inData );

if wid ~= 1 && hei ~= 1
    [Rall] = genSparseGroupingMatrix( hei, wid, 1 );
    Rblk = cell( InputData.BlkDS.blkNum, 1 );
    for i = 1:InputData.BlkDS.blkNum
        Rb = Rall(:, InputData.BlkDS.B2GMap{i});
        ins = sum(abs(Rb), 2);
        Rb(ins<2,:) = [];
        Rblk{i} = Rb;
    end
    InputData.Rblk = Rblk;
else
    InputData.Rblk = [];
end
InputData.scaleFactor = 1;
% add 1 to the zero entries in the data area
idx = find( InputData.BlkDS.indMap == 1 );
[I, J] = find(InputData.dataCube(:,idx)==0);
for i = 1:length(I)
    InputData.dataCube(I(i), idx(J(i)) ) = 1;
end
InputData.logFY = logfactorial_e( inData, 1e7);

if leaveOutPercent == 0
    InputData.aMatrix = InputData.BlkDS.indMap;
else
    InputData.aMatrix = ...
        leaveoutCBPatternsData_v2( InputData.BlkDS, leaveOutPercent );
end


end

