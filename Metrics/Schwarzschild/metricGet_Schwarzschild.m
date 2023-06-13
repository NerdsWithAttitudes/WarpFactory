function [metric] = metricGet_Schwarzschild(gridSize,worldCenter,rs,gridScaling)
%% METRICGET_SCHWARZSCHILD: Builds the Schwarzschild metric
%
%   INPUTS: 
%   gridSize - 1x4 array. world size in [t, x, y, z], double type.
%   
%   worldCenter - 1x4 array. world center location in [t, x, y, z], double type.
% 
%   rs - Schwarzschild radius
% 
%   gridScale - scaling of the grid in [t, x, y, z]. double type.
%
%   OUTPUTS: 
%   metric - metric struct object. 

%%

% Handle default input arguments
if nargin < 4
    gridScaling = [1,1,1,1];
end

% Check if gridSize in time is 1 and return warning if not
if gridSize(1) > 1
    error('The time grid is greater than 1, only a size of 1 can be used for the Schwarzschild solution');
end

% Assign parameters to metric struct
metric.params.gridSize = gridSize;
metric.params.worldCenter = worldCenter;
metric.params.rs = rs;

% Assign quantities to metric struct
metric.type = "metric";
metric.frame = "comoving";
metric.name = "Schwarzschild";
metric.scaling = gridScaling;
metric.coords = "cartesian";
metric.index = "covariant";
metric.date = date;


% Set Minkowski terms
metric.tensor = setMinkowski(gridSize);

% Add very small offset to mitigate divide by zero errors
epsilon = 0.0000000001;
t = 1; % Only 1 time slice
for i = 1:gridSize(2)
    for j = 1:gridSize(3)
        for k = 1:gridSize(4)

            x = i*gridScaling(2)-worldCenter(2);
            y = j*gridScaling(3)-worldCenter(3);
            z = k*gridScaling(4)-worldCenter(4);

            r = sqrt(x^2+y^2+z^2)+epsilon;
            
            % Diagonal terms
            metric.tensor{1,1}(t,i,j,k) = -(1-rs/r);
            metric.tensor{2,2}(t,i,j,k) = (x^2/(1-rs/r)+y^2 + z^2)/r^2;
            metric.tensor{3,3}(t,i,j,k) = (x^2 + y^2/(1-rs/r) + z^2)/r^2;
            metric.tensor{4,4}(t,i,j,k) = (x^2 + y^2 + z^2/(1-rs/r))/r^2;
            
            % dxdy cross terms
            metric.tensor{2,3}(t,i,j,k) = rs/(r^3-r^2*rs)*x*y;
            metric.tensor{3,2}(t,i,j,k) = metric.tensor{2,3}(1,i,j,k);
            
            % dxdz cross terms
            metric.tensor{2,4}(t,i,j,k) = rs/(r^3-r^2*rs)*x*z;
            metric.tensor{4,2}(t,i,j,k) = metric.tensor{2,4}(1,i,j,k);
            
            % dydz cross terms
            metric.tensor{3,4}(t,i,j,k) = rs/(r^3-r^2*rs)*y*z;
            metric.tensor{4,3}(t,i,j,k) = metric.tensor{3,4}(1,i,j,k);
        end
    end
end


