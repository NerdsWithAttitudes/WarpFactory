function [transformedEnergyTensor] = doFrameTransfer(metric, energyTensor, frame, tryGPU)

%% GETFRAMETRANSFER: Transforms the energy tensor into selected frames
%
%   INPUTS:
%   metric - Metric struct
%
%   energyTensor - Energy struct
%
%   frame - Frame to transform the tensor to. Only 'Eulerian' is supported.
%
%   tryGPU - A flag on whether or not to use GPU computation (0=no, 1=yes)
%
%
%   OUTPUTS:
%   transformedEnergyTensor - Transformed energy struct

%%

% Handle default input arguments
if nargin < 4
    tryGPU = 0;
end
transformedEnergyTensor = energyTensor;
transformedEnergyTensor.tensor = cell(4, 4);

% Check metric and tensor format is correct
if ~verifyTensor(metric, 1)
    error("Metric is not verified. Please verify metric using verifyTensor(metric).")
end
if ~verifyTensor(energyTensor, 1)
    error("Stress-energy is not verified. Please verify Stress-energy tensor using verifyTensor(energyTensor).")
end

if strcmpi(frame, "Eulerian") && ~(isfield(energyTensor, 'frame') && strcmpi(energyTensor.frame, 'Eulerian'))
    % Convert to covariant (lower) index
    energyTensor = changeTensorIndex(energyTensor, "covariant", metric);

    % Convert from cell to array
    arrayEnergyTensor = tensorCell2Array(energyTensor, tryGPU);
    arrayMetricTensor = tensorCell2Array(metric, tryGPU);

    % Do transformations at each point in space
    M = getEulerianTransformationMatrix(arrayMetricTensor, metric.coords);
    M = permute(M, [5, 6, 1, 2, 3, 4]);
    arrayEnergyTensor = permute(arrayEnergyTensor, [5, 6, 1, 2, 3, 4]);

    transformedTempTensor.tensor = pagemtimes(pagemtimes(M, 'transpose', arrayEnergyTensor, 'none'), M);

    % Convert array tensor into cell format
    z = size(transformedTempTensor.tensor);
    for i = 1:4
        for j = 1:4
            transformedEnergyTensor.tensor{i, j} = reshape(transformedTempTensor.tensor(i, j, :, :, :, :), [z(3:end), 1]);
        end
    end

    % Tranform to contravariant T^{0, i} = -T_{0, i}
    for i = 2:4
        transformedEnergyTensor.tensor{1, i} = -transformedEnergyTensor.tensor{1, i};
        transformedEnergyTensor.tensor{i, 1} = -transformedEnergyTensor.tensor{i, 1};
    end

    % Update the tensor metadata
    transformedEnergyTensor.frame = "Eulerian";
    transformedEnergyTensor.index = "contravariant";

else
    warning("Frame not found")
end

end
