function [metric] = metricGet_ModifiedTimeComoving(gridSize,worldCenter,v, R, sigma, A, gridScaling)
%% METRICGET_MODIFIEDTIMECOMOVING: Builds the Modified Time metric in Galilean comoving frame
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
%   A - lapse rate modification, double type.
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
    error('The time grid is greater than 1, only a size of 1 can be used in comoving');
end

% Assign parameters to metric struct
metric.params.gridSize = gridSize;
metric.params.worldCenter = worldCenter;
metric.params.velocity = v;
metric.params.R = R;
metric.params.sigma = sigma;
metric.params.A = A;

% Assign parameters to metric struct
metric.type = "metric";
metric.name = "Modified Time Comoving";
metric.scaling = gridScaling;
metric.coords = "cartesian";
metric.index = "covariant";
metric.date = date;


% Declare a Minkowski space
metric.tensor = setMinkowski(gridSize);

% Add the Modified Time changes
t = 1; % only one timeslice is used
for i = 1:gridSize(2)
    for j = 1:gridSize(3)
        for k = 1:gridSize(4)

            x = i*gridScaling(1)-worldCenter(1);
            y = j*gridScaling(2)-worldCenter(2);
            z = k*gridScaling(3)-worldCenter(3);

            % Find the radius from the center of the bubble
            r = (x^2 + y^2 + z^2)^(1/2);

            % Find shape function at this point in r
            fs = shapeFunction_Alcubierre(r,R,sigma);

            % Add alcubierre term to dxdt
            metric.tensor{1,2}(t,i,j,k) = v*(1-fs);
            metric.tensor{2,1}(t,i,j,k) = metric.tensor{1,2}(t,i,k,k);

            % Add dt term modification
            metric.tensor{1,1}(t,i,j,k) = -((1-fs)+fs/A)^2+(fs*v)^2;        
        end
    end
end



