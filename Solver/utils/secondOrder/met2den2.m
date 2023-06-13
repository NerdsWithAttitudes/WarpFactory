function [energyDensity] = met2den2(metricTensor,delta,units)
%MET2DEN Coverts a catesian metric tensor to the corresponding energy 
%   density tensor using the Einstien Field Equations
%   Takes an input of a 4x4 cell array as the metric tensor and outputs a
%   4x4 cell array as the energy density tensor
%
%   INPUT: 4x4 cell array. Elements of cell array are 4-D matricies of
%   double type
%
%
%   OUPUT: 4x4 cell array. Elements of cell array are 4-D matricies of
%   double type

switch nargin
    case 1
        delta = [1,1,1,1];
        units = [1,1,1];
    case 2
        units = [1,1,1];
end

% Metric tensor and its inverse
gl = metricTensor;
gu = c4Inv2(gl);

% Calculate the Ricci tensor
R_munu = ricciT2(gu,gl,delta);

% Calculate the Ricci scalar
R = ricciS2(R_munu,gu);

% Calculate Einstien tensor
E = einT2(R_munu,R,gl);

% Calculate Energy density
energyDensity = einE2(E,gu,units);

end


