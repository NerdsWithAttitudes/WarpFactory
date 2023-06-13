function [energyDensity] = met2den(gl, delta)

%% MET2DEN: Coverts a catesian metric tensor to the corresponding energy density tensor using the Einstien Field Equations
%
%   INPUTS:
%   gl - 4x4 cell array. Elements of cell array are 4-D matricies of
%   double type
%
%   delta - 1x4 array of the uniform delta step size in each coordinate
%   direction
%
%
%   OUTPUTS:
%   energyDensity - 4x4 cell array. Elements of cell array are 4-D matricies of
%   double type

%%

% Handle default input arguments
if nargin < 2
    delta = [1, 1, 1, 1];
end

% Calcualte metric tensor inverse
gu = c4Inv(gl);

% Calculate the Ricci tensor
R_munu = ricciT(gu, gl, delta);

% Calculate the Ricci scalar
R = ricciS(R_munu, gu);

% Calculate Einstien tensor
E = einT(R_munu, R, gl);

% Calculate Energy density
energyDensity = einE(E, gu);

end
