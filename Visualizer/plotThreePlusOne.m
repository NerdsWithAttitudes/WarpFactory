function plotThreePlusOne(metric, slicedPlanes, sliceLocations, alpha)

%% PLOTTHREEPLUSONE: Plots the 3+1 elements of the METRIC tensor based the slice plane
%
%   INPUTS:
%   metric - metric tensor struct object.
%
%   slicedPlanes - Coordinates that are sliced [coords1, coords2], index values from 1 to 4, double type.
%   If you want the resulting slice to be in the X-Y plane for example, input [1, 4]
%
%   sliceLocation - Location of the slice in [coords1, coords2], double type.
%
%   alpha - Alpha value of the surface grid display from 0 to 1, double type.
%
%   OUTPUTS:
%   [none]

%%

% Handle input arguments
if nargin < 2
    slicedPlanes = [1, 4]; % Assume X,Y plane
end
if nargin < 3
    s = size(metric.tensor{1, 1});
    sliceLocations = round((s+1)./2); % Assume center
end
if nargin < 4
    alpha = 0.2;
end

% Check that tensor is a metric
if strcmpi(metric.type, "metric")
    error("Must provide a metric object.");
end

% Verify tensor
if ~verifyTensor(metric, 1)
    error("Metric is not verified. Please verify metric using verifyTensor(metric).")
end

% Check that the sliced planes are different
if slicedPlanes(1) == slicedPlanes(2)
    error("Selected planes must not be the same, select two different planes to slice along.")
end

% Round sliceLocations
sliceLocations = round(sliceLocations);

% Check that the sliceLocations are inside the world
if sliceLocations(1) < 1 || sliceLocations(2) < 1 || sliceLocations(1) > size(tensor.tensor, slicedPlanes(1)) || sliceLocations(2) > size(tensor.tensor, slicedPlanes(2))
    error('sliceLocations are outside the world.')
end

% Check that the coords are cartesian
if strcmpi(metric.coords, "cartesian")
    [alpha_lapse, betaDown, gammaDown, ~, ~] = ThreePlusOneDecomposer(metric);

    [xLabelText, yLabelText] = labelCartesianAxis(slicedPlanes);
    idx = getSliceData(slicedPlanes, sliceLocations, metric);

    % Plot alpha
    titleText = "\alpha";
    plotComponent(squeeze(alpha_lapse(idx{1}, idx{2}, idx{3}, idx{4}))', titleText, xLabelText, yLabelText, alpha)

    % Plot beta
    for i = 1:3
        titleText = "\beta_"+num2str(i);
        plotComponent(squeeze(betaDown{i}(idx{1}, idx{2}, idx{3}, idx{4}))', titleText, xLabelText, yLabelText, alpha)
    end

    % Plot gamma
    c = [1, 1; 1, 2; 1, 3; 2, 2; 2, 3; 3, 3];

    for i = 1:size(c, 1)
        titleText = "\gamma_{"+num2str(c(i, 1))+num2str(c(i, 2))+"}";
        plotComponent(squeeze(gammaDown{c(i, 1), c(i, 2)}(idx{1}, idx{2}, idx{3}, idx{4}))', titleText, xLabelText, yLabelText, alpha)
    end

else
    error('Unknown coordinate system, must be: "cartesian"')
end


end
