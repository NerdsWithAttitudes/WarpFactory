function metric = setMinkowski(gridSize)
%% SETMINKOWSKI: Builds metric terms for a flat Minkowski space
%
%   INPUTS:
%   gridSize - World size in [t,x,y,z]
%
%   metric - Metric struct
%
%   OUTPUTS:
%   tensor - The metric tensor as a 4x4 cell of 4D arrays.

%%

% dt^2 term
metric{1,1} = -ones(gridSize);


% Non-time diagonal terms
metric{2,2} = ones(gridSize);
metric{3,3} = ones(gridSize);
metric{4,4} = ones(gridSize);

% Cross terms
metric{1,2} = zeros(gridSize);
metric{2,1} = zeros(gridSize);
metric{1,3} = zeros(gridSize);
metric{3,1} = zeros(gridSize);
metric{2,3} = zeros(gridSize);
metric{3,2} = zeros(gridSize);
metric{2,4} = zeros(gridSize);
metric{4,2} = zeros(gridSize);
metric{3,4} = zeros(gridSize);
metric{4,3} = zeros(gridSize);
metric{1,4} = zeros(gridSize);
metric{4,1} = zeros(gridSize);

   


end