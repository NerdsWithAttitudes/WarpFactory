function [alpha, beta, gamma] = setMinkowskiThreePlusOne(gridSize)
%% SETMINKOWSKI: Returns the 3+1 format for flat space
%
%   INPUTS:
%   gridSize - World size in [t,x,y,z]
%
%   OUTPUTS:
%   alpha - Lapse rate 4D array
%
%   beta - Shift vector, 1x3 cell of 4D arrays
%
%   gamma - Spatial terms, 3x3 cell of 4D arrays.

%%

alpha = ones(gridSize);

beta = cell(1,3);
for i = 1:3
    beta{i} = zeros(gridSize);
end

gamma = cell(3,3);
for i = 1:3
    for j = 1:3
        if i == j
            gamma{i,j} = ones(gridSize);
        else
            gamma{i,j} = zeros(gridSize);
        end
    end
end
 
end