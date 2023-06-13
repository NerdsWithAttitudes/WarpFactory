function energy = getEnergyTensor(metric, tryGPU, diffOrder)

%% GETENERGYTENSOR: Converts the metric into the stress energy tensor
%
%   INPUTS:
%   metric - A metric struct
%
%   tryGPU - A flag on whether or not to use GPU computation (0=no, 1=yes)
%
%   diffOrder - Order of finite difference, either 'second' or 'fourth'
%
%   OUTPUTS:
%   energy - energy tensor struct

%%

% Handle default input arguments
if nargin < 2
    tryGPU = 0;
end
if nargin < 3
    diffOrder = 'fourth';
end

% Check that the metric is verified and covariant
if ~verifyTensor(metric, 1)
    error("Metric is not verified. Please verify metric using verifyTensor(metric).")
end
if ~strcmpi(metric.index, "covariant")
    metric = changeTensorIndex(metric, "covariant");
    fprintf("Changed metric from %s index to %s index\n", metric.index, "covariant")
end

% Use GPU for computation
if tryGPU
    % Convert metric to GPU Array
    metricTensorGPU = cell(4, 4);
    for i = 1:4
        for j = 1:4
            metricTensorGPU{i, j} = gpuArray(metric.tensor{i, j});
        end
    end

    % Compute on GPU
    metric.scaling = gpuArray(metric.scaling);
    if strcmp(diffOrder, 'fourth')
        enDenGPU = met2den(metricTensorGPU, metric.scaling);
    elseif strcmp(diffOrder, 'second')
        enDenGPU = met2den2(metricTensorGPU, metric.scaling);
    else
        error("Order Flag Not Specified Correctly. Options: 'fourth' or 'second'")
    end

    % Gather results from GPU
    energyTensor = cell(4, 4);
    for i = 1:4
        for j = 1:4
            energyTensor{i, j} = gather(enDenGPU{i, j});
        end
    end

else
    % Compute on CPU
    if strcmp(diffOrder, 'fourth')
        energyTensor = met2den(metric.tensor, metric.scaling);
    elseif strcmp(diffOrder, 'second')
        energyTensor = met2den2(metric.tensor, metric.scaling);
    else
        error("Order Flag Not Specified Correctly. Options: 'fourth' or 'second'")
    end
end

% Assign struct values
energy.type = "Stress-Energy";
energy.tensor = energyTensor;
energy.coords = metric.coords;
energy.index = "contravariant";
energy.order = diffOrder;
energy.name = metric.name;
energy.date = date;

end
