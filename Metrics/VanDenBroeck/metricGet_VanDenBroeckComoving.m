function [metric] = metricGet_VanDenBroeckComoving(gridSize, worldCenter, v, R1, sigma1, R2, sigma2, A, gridScale)
%% METRICGET_VANDENBROECKCOMOVING: Builds the Van Den Broeck metric in Galilean comoving frame
%
%   INPUTS: 
%   gridSize - 1x4 array. world size in [t, x, y, z], double type.
%   
%   worldCenter - 1x4 array. world center location in [t, x, y, z], double type.
% 
%   v - speed of the warp drive in factors of c, along the x direction, double type.
% 
%   R1 - spatial expansion radius of the warp bubble, double type.
% 
%   sigma1 - width factor of the spatial expansion transition
%
%   R2 - shift vector radius of the warp bubble, double type.
%
%   sigma2 - width factor of the shift vector transition
%
%   A - spatial expansion factor, double type.
% 
%   gridScale - scaling of the grid in [t, x, y, z]. double type.
% 
%   OUTPUTS: 
%   metric - metric struct object. 

%%

% Handle default input arguments
if nargin < 9
    gridScale = [1,1,1,1];
end

% Check if gridSize in time is 1 and return warning if not
if gridSize(1) > 1
    error('The time grid is greater than 1, only a size of 1 can be used in comoving');
end


% Assign parameters to metric struct
metric.params.gridSize = gridSize;
metric.params.worldCenter = worldCenter;
metric.params.velocity = v*(1+A)^2;
metric.params.R1 = R1;
metric.params.sigma1 = sigma1;
metric.params.R2 = R2;
metric.params.sigma2 = sigma2;
metric.params.A = A;

% Assign quantities to metric struct
metric.type = "metric";
metric.name = 'Van Den Broeck Comoving';
metric.scaling = gridScale;
metric.coords = 'cartesian';
metric.index = "covariant";
metric.date = date;

%% Build Metric

% Declare a Minkowski space
metric.tensor = setMinkowski(gridSize);

% Van Den Brock modification
t = 1; % only one timeslice is used
for i = 1:gridSize(2)
    for j = 1:gridSize(3)
        for k = 1:gridSize(4)

            x = i*gridScale(2)-worldCenter(2);
            y = j*gridScale(3)-worldCenter(3);
            z = k*gridScale(4)-worldCenter(4);
            
            % Find the radius from the center of the bubble
            r = sqrt(x^2 + y^2 + z^2);

            % Define the B function value in Van Den Broeck
            B = 1+shapeFunction_Alcubierre(r,R1,sigma1)*A;

            % Define the f function value in Van Den Broeck
            fs = shapeFunction_Alcubierre(r,R2,sigma2)*v;
    
            % Assign fs and B to the proper terms
            metric.tensor{2,2}(t,i,j,k) = B^2;
            metric.tensor{3,3}(t,i,j,k) = B^2;
            metric.tensor{4,4}(t,i,j,k) = B^2;

            metric.tensor{1,2}(t,i,j,k) = B^2*(v-fs);
            metric.tensor{2,1}(t,i,j,k) = metric.tensor{1,2}(t,i,j,k);

            metric.tensor{1,1}(t,i,j,k) = -(1-B^2*fs^2);

        end
    end
end


end