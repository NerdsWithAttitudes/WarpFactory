function [paths] = getMomentumFlowLines(energyTensor, startPoints, stepSize, maxSteps, scaleFactor)

%% GETMOMENTUMFLOWLINES: Gets the momentum flow lines for an energy tensor
%
%   INPUTS:
%   energyTensor - Energy struct
%
%   startPoints - 1x3 cell array of the start points of flowlines
%   startPoints{1} = X;
%   startPoints{2} = Y;
%   startPoints{3} = Z;
%
%   stepSize - Step size of the flowline propagation
%
%   maxSteps - The max number of propagation steps to run
%
%   scaleFactor - The scaling factor that multiplies the momentum density
%
%
%   OUTPUTS:
%   paths - 1xN cell array containing N paths. The path in each cell is an Mx3 array.
%

%%

% Check that the energyTensor is contravariant
if ~strcmpi(energyTensor.index, "contravariant")
    error('Energy tensor for momentum flowlines should be contravariant.')
end

% Load in the momentum data
Xmom = squeeze(energyTensor.tensor{1, 2}) * scaleFactor;
Ymom = squeeze(energyTensor.tensor{1, 3}) * scaleFactor;
Zmom = squeeze(energyTensor.tensor{1, 4}) * scaleFactor;

% Reshape the starting points X, Y, and Z
StrPtsX = reshape(startPoints{1}, 1, numel(startPoints{1}));
StrPtsY = reshape(startPoints{2}, 1, numel(startPoints{2}));
StrPtsZ = reshape(startPoints{3}, 1, numel(startPoints{3}));

% Make the paths
paths = cell(1, length(StrPtsX));
for j = 1:length(StrPtsX)
    Pos = zeros(maxSteps, 3);
    Pos(1, :) = [StrPtsX(j), StrPtsY(j), StrPtsZ(j)];

    for i = 1:maxSteps
        % Check if the particle is outside the world
        if sum(isnan(Pos(i, :))) > 0 || (floor(Pos(i, 1)) <= 1 || ceil(Pos(i, 1)) >= size(Xmom, 1)) || ...
                (floor(Pos(i, 2)) <= 1 || ceil(Pos(i, 2)) >= size(Xmom, 2)) || ...
                (floor(Pos(i, 3)) <= 1 || ceil(Pos(i, 3)) >= size(Xmom, 3))
            break;
        end

        % Interpolate the momentum
        xMomentum = trilinearInterp(Xmom, Pos(i, :));
        yMomentum = trilinearInterp(Ymom, Pos(i, :));
        zMomentum = trilinearInterp(Zmom, Pos(i, :));

        % Propagate position
        Pos(i+1, 1) = Pos(i, 1)+(xMomentum) * stepSize;
        Pos(i+1, 2) = Pos(i, 2)+(yMomentum) * stepSize;
        Pos(i+1, 3) = Pos(i, 3)+(zMomentum) * stepSize;
    end

    paths{j} = Pos(1:i-1, :);
end

end