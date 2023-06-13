function metric = metricGet_Minkowski(gridSize,gridScaling)
%% METRICGET_MINKOWSKI: Builds a Minkowski metric
%
%   INPUTS: 
%   gridSize - 1x4 array. world size in [t, x, y, z], double type. 
% 
%   gridScale - scaling of the grid in [t, x, y, z]. double type.
% 
%   OUTPUTS: 
%   metric - metric struct object. 

%%


% Handle default input arguments
if nargin < 2
    gridScaling = [1,1,1,1];
end

% Assign quantities to metric struct
metric.type = "metric";
metric.name = "Minkowski";
metric.scaling = gridScaling;
metric.coords = "cartesian";
metric.index = "covariant";
metric.date = date;


% dt^2 term
metric.tensor{1,1} = -ones(gridSize);

% Non-time diagonal terms
metric.tensor{2,2} = ones(gridSize);
metric.tensor{3,3} = ones(gridSize);
metric.tensor{4,4} = ones(gridSize);

% Cross terms
metric.tensor{1,2} = zeros(gridSize);
metric.tensor{2,1} = zeros(gridSize);
metric.tensor{1,3} = zeros(gridSize);
metric.tensor{3,1} = zeros(gridSize);
metric.tensor{2,3} = zeros(gridSize);
metric.tensor{3,2} = zeros(gridSize);
metric.tensor{2,4} = zeros(gridSize);
metric.tensor{4,2} = zeros(gridSize);
metric.tensor{3,4} = zeros(gridSize);
metric.tensor{4,3} = zeros(gridSize);
metric.tensor{1,4} = zeros(gridSize);
metric.tensor{4,1} = zeros(gridSize);


end