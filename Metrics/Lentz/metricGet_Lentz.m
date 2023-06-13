function [metric] = metricGet_Lentz(gridSize, worldCenter, v, scale, gridScale)
%% METRICGET_LENTZ: Builds the Lentz metric
%
%   INPUTS: 
%   gridSize - 1x4 array. world size in [t, x, y, z], double type.
%   
%   worldCenter - 1x4 array. world center location in [t, x, y, z], double type.
% 
%   v - speed of the warp drive in factors of c, along the x direction, double type.
%
%   scale - the sizing factor of the Lentz soliton template
% 
%   gridScale - scaling of the grid in [t, x, y, z]. double type.
%
%   OUTPUTS: 
%   metric - metric struct object. 

%%

% Handle default input arguments
if nargin < 4
    scale = max(gridSize(2:4))/7;
end
if nargin < 5
    gridScale = [1,1,1,1];
end

% Assign parameters to metric struct
metric.params.gridSize = gridSize;
metric.params.worldCenter = worldCenter;
metric.params.velocity = v;

% Assign quantities to metric struct
metric.type = "metric";
metric.name = "Lentz";
metric.scaling = gridScale;
metric.coords = "cartesian";
metric.index = "covariant";
metric.date = date;


% Declare a Minkowski space
[alpha, beta, gamma] = setMinkowskiThreePlusOne(gridSize);

% Lentz Soliton Terms
for i = 1:gridSize(2)
    for j = 1:gridSize(3)
        for k = 1:gridSize(4)
            
            x = i*gridScale(2)-worldCenter(2);
            y = j*gridScale(3)-worldCenter(3);

            for t = 1:gridSize(1)
                % Determine the x offset of the center of the bubble,
                % centered in time
                xs = (t*gridScale(1)-worldCenter(1))*v*c;

                xp = x-xs;

                % Get Lentz template values
                [WFX, WFY] = getWarpFactorByRegion(xp,y,scale);
        
                % Assign dxdt term
                beta{1}(t,i,j,k) = -WFX*v;

                % Assign dydt term
                beta{2}(t,i,j,k) = WFY*v;
            end
        end
    end
end

% Make tensor from the 3+1 functions
metric.tensor = threePlusOneBuilder(alpha,beta,gamma);

end



function [WFX, WFY] = getWarpFactorByRegion(xIn,yIn,sizeScale)
x = xIn;
y = abs(yIn);
WFX = 0;
WFY = 0;

% Lentz shift vector template
if (x >= sizeScale && x <= 2*sizeScale) && (x-sizeScale >= y)
    WFX = -2;
    WFY = 0;
elseif (x > sizeScale && x <= 2*sizeScale) && (x-sizeScale <= y) && (-y+3*sizeScale >= x)
    WFX = -1;
    WFY = 1;
elseif (x > 0 && x <= sizeScale) && (x+sizeScale > y) && (-y+sizeScale < x)
    WFX = 0;
    WFY = 1;
elseif (x > 0 && x <= sizeScale) && (x+sizeScale <= y) && (-y+3*sizeScale >= x)
    WFX = -0.5;
    WFY = 0.5;
elseif (x > -sizeScale && x <= 0) && (-x+sizeScale < y) && (-y+3*sizeScale >= -x)
    WFX = 0.5;
    WFY = 0.5;
elseif (x > -sizeScale && x <= 0) && (x+sizeScale <= y) && (-y+sizeScale >= x)
    WFX = 1;
    WFY = 0;
elseif (x >= -sizeScale && x <= sizeScale) && (x+sizeScale > y)
    WFX = 1;
    WFY = 0;
end

WFY = sign(yIn)*WFY;

end



