function plotTensor(tensor, alpha, slicedPlanes, sliceLocations)

%% PLOTTENSOR: Plots the unique elements of the tensor based the slice plane
%
%   INPUTS:
%   tensor - Tensor struct object either metric of stress-energy.
%
%   alpha - Alpha value of the surface grid display from 0 to 1, double type.
%
%   slicedPlanes - Coordinates that are sliced [coords1, coords2], index values from 1 to 4, double type.
%   If you want the resulting slice to be in the X-Y plane for example, input [1, 4]
%
%   sliceLocation - Location of the slice in [coords1, coords2], double type.
%
%   OUTPUTS:
%   [none]

%%

% Handle default input arguments
if nargin < 2
    alpha = 0.2;
end
if nargin < 3
    slicedPlanes = [1, 4]; % Assume X,Y plane
end
if nargin < 4
    s = size(tensor.tensor{1, 1});
    sliceCenters = round((s+1)./2); % Assume center
    sliceLocations(1) = sliceCenters(slicedPlanes(1));
    sliceLocations(2) = sliceCenters(slicedPlanes(2));
end

% Verify tensor
if ~verifyTensor(tensor, 1)
    error("Tensor is not verified. Please verify tensor using verifyTensor(tensor).")
end

% Check that the sliced planes are different
if slicedPlanes(1) == slicedPlanes(2)
    error("Selected planes must not be the same, select two different planes to slice along.")
end

% Round sliceLocations
sliceLocations = round(sliceLocations);

% Check that the sliceLocations are inside the world
if sliceLocations(1) < 1 || sliceLocations(2) < 1 || sliceLocations(1) > size(tensor.tensor{1, 1}, slicedPlanes(1)) || sliceLocations(2) > size(tensor.tensor{1, 1}, slicedPlanes(2))
    sliceLocations(1)
    sliceLocations(2)   
    size(tensor.tensor{1, 1}, slicedPlanes(1)) 
    size(tensor.tensor{1, 1}, slicedPlanes(2))
    error('sliceLocations are outside the world.')
end

% Check tensor type
if strcmpi(tensor.type, "Metric")
    titleCharacter = "g";
elseif strcmpi(tensor.type, "Stress-Energy")
    titleCharacter = "T";
end

% Check tensor index
if strcmpi(tensor.index, "covariant")
    titleAugment1 = "_{";
    titleAugment2 = "";
elseif strcmpi(tensor.index, "contravariant")
    titleAugment1 = "^{";
    titleAugment2 = "";
elseif strcmpi(tensor.index, "mixedupdown")
    titleAugment1 = "^{";
    titleAugment2 = "}_{ ";
elseif strcmpi(tensor.index, "mixeddownup")
    titleAugment1 = "_{";
    titleAugment2 = "}^{ ";
end

% Check that the coords are cartesian
if strcmpi(tensor.coords, "cartesian")
    [xLabelText, yLabelText] = labelCartesianAxis(slicedPlanes);

    if strcmpi(tensor.index, "mixedupdown") || strcmpi(tensor.index, "mixeddownup")
        c1 = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4];
        c2 = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4];
    else
        c1 = [1, 1, 1, 1, 2, 3, 4, 2, 2, 3];
        c2 = [1, 2, 3, 4, 2, 3, 4, 3, 4, 4];
    end

    idx = getSliceData(slicedPlanes, sliceLocations, tensor);

    for i = 1:length(c1)
        plotComponent(squeeze(tensor.tensor{c1(i), c2(i)}(idx{1}, idx{2}, idx{3}, idx{4}))', titleCharacter+titleAugment1+num2str(c1(i))+titleAugment2+num2str(c2(i))+"}", xLabelText, yLabelText, alpha)
    end
else
    error('Unknown coordinate system, must be: "cartesian"')
end


end