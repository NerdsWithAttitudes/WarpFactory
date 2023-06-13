function [metric] = metricGet_Alcubierre(gridSize,worldCenter,v,R,sigma,gridScale)
%% METRICGET_ALCUBIERRE: Builds the Alcubierre metric
%
%   INPUTS: 
%   gridSize - 1x4 array. world size in [t, x, y, z], double type.
%   
%   worldCenter - 1x4 array. world center location in [t, x, y, z], double type.
% 
%   v - speed of the warp drive in factors of c, along the x direction, double type.
% 
%   R - radius of the warp bubble, double type.
% 
%   sigma - thickness parameter of the bubble, double type.
% 
%   gridScale - scaling of the grid in [t, x, y, z]. double type.
%
%   OUTPUTS: 
%   metric - metric struct object. 

%%


% Handle default input arguments
if nargin < 6
    gridScale = [1,1,1,1];
end

% Assign parameters to metric struct
metric.params.gridSize = gridSize;
metric.params.worldCenter = worldCenter;
metric.params.velocity = v;
metric.params.R = R;
metric.params.sigma = sigma;

% Assign quantities to metric struct
metric.type = "metric";
metric.name = 'Alcubierre';
metric.scaling = gridScale;
metric.coords = "cartesian";
metric.index = "covariant";
metric.date = date;

% Declare a Minkowski space
[alpha, beta, gamma] = setMinkowskiThreePlusOne(gridSize);

% Add the Alcubierre modification
for i = 1:gridSize(2)
    for j = 1:gridSize(3)
        for k = 1:gridSize(4)
        
            % Find grid center x, y, z
            x = i*gridScale(2)-worldCenter(2);
            y = j*gridScale(3)-worldCenter(3);
            z = k*gridScale(4)-worldCenter(4);

            for t = 1:gridSize(1)
                % Determine the x offset of the center of the bubble,
                % centered in time
                xs = (t*gridScale(1)-worldCenter(1))*v*c;

                % Find the radius from the center of the bubble
                r = ((x-xs)^2 + y^2 + z^2)^(1/2);

                % Find shape function at this point in r
                fs = shapeFunction_Alcubierre(r,R,sigma);

                % Add alcubierre modification to shift vector along x
                beta{1}(t,i,j,k) = -v*fs;
            end
        end
    end
end

% Make tensor from the 3+1 functions
metric.tensor = threePlusOneBuilder(alpha,beta,gamma);
end