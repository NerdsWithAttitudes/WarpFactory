function [map, vec, vectorFieldOut] = getEnergyConditions(energyTensor, metric, condition, numAngularVec, numTimeVec, returnVec, tryGPU)

%% GETENERGYCONDITIONS: Function to get the energy conditions of an energy tensor
%
%   INPUTS:
%   energyTensor - Energy struct
%
%   metric - Metric struct
%
%   condition - What energy condition to evaluate. Either "Null", "Weak", "Strong", or "Dominant"
%
%   numAngularVec - Number of equally spaced spatial vectors to evaluate
%
%   numTimeVec - Number of equally spaced temporal shells to evaluate
%
%   returnVec - A flag on whether or not to return all evaluations and their vectors (0=no, 1=yes)
%
%   tryGPU - A flag on whether or not to use GPU computation (0=no, 1=yes)
%
%
%   OUTPUTS:
%   map - The most violating evaluation at every point in the spacetime
%
%   vec - The evaluations for every vector at every point in spacetime
%   (only returned when returnVec is 1)
%
%   vectorFieldOut - The vector field used at each point in the spacetime
%   to evaluate the enegy conditions (only returned when returnVec is 1)

%%

% Handle default input arguments
if nargin < 4
    numAngularVec = 100;
end
if nargin < 5
    numTimeVec = 10;
end
if nargin < 6
    returnVec = 0;
end
if nargin < 7
    tryGPU = 0;
end

% Check if correct conditions input
if ~(strcmpi(condition, "Null") || strcmpi(condition, "Weak") || strcmpi(condition, "Dominant") || strcmpi(condition, "Strong"))
    error('Incorrect energy condition input, use either: "Null", "Weak", "Dominant", "Strong"')
end

% Return warning for any coordinate system not cartesian
if ~strcmpi(metric.coords, 'cartesian')
    warning('Evaluation not verified for coordinate systems other than Cartesian!')
end

% Check tensor formats are correct
if ~verifyTensor(metric, 1)
    error("Metric is not verified. Please verify metric using verifyTensor(metric).")
end
if ~verifyTensor(energyTensor, 1)
    error("Stress-energy is not verified. Please verify stress-energy using verifyTensor(EnergyTensor).")
end

% Convert arrays into GPU if needed
if tryGPU
    energyTensorGPU = energyTensor;
    MetricGPU = metric;
    for i = 1:4
        for j = 1:4
            energyTensorGPU.tensor{i, j} = gpuArray(energyTensor.tensor{i, j});
            MetricGPU.tensor{i, j} = gpuArray(metric.tensor{i, j});
        end
    end
    energyTensor = energyTensorGPU;
    metric = MetricGPU;
end

% Get size of the spacetime
[a, b, c, d] = size(metric.tensor{1, 1});

% Convert energy tensor into the local inertial frame if not eulerian
energyTensor = doFrameTransfer(metric, energyTensor, "Eulerian", tryGPU);

%% Build vector fields
if strcmpi(condition, "Null") || strcmpi(condition, "Dominant")
    type = "nulllike";
elseif strcmp(condition, "Weak") || strcmp(condition, "Strong")
    type = "timelike";
end
vecField = generateUniformField(type, numAngularVec, numTimeVec, tryGPU);

% Declare variables to be determined in eval of energy conditions
if isgpuarray(metric.tensor{1, 1})
    vecField = gpuArray(vecField);
    map = nan(a, b, c, d, 'gpuArray');
    if returnVec == 1
        vec = zeros(a, b, c, d, numAngularVec, numTimeVec, 'gpuArray');
    end
else
    map = nan(a, b, c, d);
    if returnVec == 1
        vec = zeros(a, b, c, d, numAngularVec, numTimeVec);
    end
end

%% Find energy conditions

% Null energy condition
if strcmpi(condition, "Null")
    energyTensor = changeTensorIndex(energyTensor, "covariant", metric); % double check that it is covariant

    for ii = 1:numAngularVec
        if isgpuarray(metric.tensor{1, 1})
            temp = zeros(a, b, c, d, 'gpuArray');
        else
            temp = zeros(a, b, c, d);
        end
        for mu = 1:4
            for nu = 1:4
                temp = temp+energyTensor.tensor{mu, nu} * vecField(mu, ii) * vecField(nu, ii);
            end
        end
        map = min(map, temp);
        if returnVec == 1
            vec(:, :, :, :, ii) = temp;
        end
    end

    % Weak energy condition
elseif strcmpi(condition, "Weak")
    energyTensor = changeTensorIndex(energyTensor, "covariant", metric); % double check that it is covariant
    for jj = 1:numTimeVec
        for ii = 1:numAngularVec
            if isgpuarray(metric.tensor{1, 1})
                temp = zeros(a, b, c, d, 'gpuArray');
            else
                temp = zeros(a, b, c, d);
            end

            for mu = 1:4
                for nu = 1:4
                    temp = temp+energyTensor.tensor{mu, nu} * vecField(mu, ii, jj) * vecField(nu, ii, jj);
                end
            end
            map = min(map, temp);
            if returnVec == 1
                vec(:, :, :, :, ii, jj) = temp;
            end
        end
    end

    % Dominant energy condition
elseif strcmpi(condition, "Dominant")
    % Build minkowski reference metric
    metricMinkowski = metricGet_Minkowski([a, b, c, d]);
    metricMinkowski = changeTensorIndex(metricMinkowski, "covariant"); %make sure it is covariant

    energyTensor = changeTensorIndex(energyTensor, "mixedupdown", metricMinkowski); % convert to mixed up down with minkowski

    for ii = 1:numAngularVec
        if isgpuarray(metric.tensor{1, 1})
            temp = zeros(a, b, c, d, 4, 'gpuArray');
        else
            temp = zeros(a, b, c, d, 4);
        end
        for mu = 1:4
            for nu = 1:4
                temp(:, :, :, :, mu) = temp(:, :, :, :, mu)-energyTensor.tensor{mu, nu} * vecField(nu, ii);
            end
        end

        % Wrap into vector struct.
        vector.field = {temp(:, :, :, :, 1), temp(:, :, :, :, 2), temp(:, :, :, :, 3), temp(:, :, :, :, 4)};
        vector.index = "contravariant";
        vector.type = "4-vector";

        % Find inner product to determine if timelike or null
        diff = getInnerProduct(vector, vector, metricMinkowski);
        diff = sign(diff) .* sqrt(abs(diff));

        map = max(map, diff);
        if returnVec == 1
            vec(:, :, :, :, ii) = diff;
        end
    end

    % Flip sign of the dominant energy condition to better align with evaluations of
    % other conditions (i.e. negative is violating)
    map = -map;
    if returnVec == 1
        vec = -vec;
    end

    % Strong energy condition
elseif strcmpi(condition, "Strong")
    % Build minkowski reference metric
    metricMinkowski = metricGet_Minkowski([a, b, c, d]);
    metricMinkowski = changeTensorIndex(metricMinkowski, "covariant");

    % Make sure the energy tensor is covariant
    energyTensor = changeTensorIndex(energyTensor, "covariant", metricMinkowski);

    % Find the trace
    ETrace = getTrace(energyTensor, metricMinkowski);

    for jj = 1:numTimeVec
        for ii = 1:numAngularVec
            if isgpuarray(metric.tensor{1, 1})
                temp = zeros(a, b, c, d, 'gpuArray');
            else
                temp = zeros(a, b, c, d);
            end

            for mu = 1:4
                for nu = 1:4
                    temp = temp+(energyTensor.tensor{mu, nu}-0.5 .* ETrace .* metricMinkowski.tensor{mu, nu}) .* vecField(mu, ii, jj) .* vecField(nu, ii, jj);
                end
            end

            map = min(map, temp);
            if returnVec == 1
                vec(:, :, :, :, ii, jj) = temp;
            end
        end
    end

else
    error('Unrecognized input energy condition, use either: "Null", "Weak", "Strong", "Dominant"')
end

% If GPU is used, gather data into arrays
if isgpuarray(metric.tensor{1, 1})
    map = gather(map);
    if returnVec == 1
        vec = gather(vec);
    end
end

% If returnVec is set, return both the vec (other places in code) and the vectorField
if returnVec == 1
    vectorFieldOut = vecField;
else
    vec = [];
    vectorFieldOut = [];
end

end
