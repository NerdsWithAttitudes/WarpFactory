function [output] = evalMetric(metric, tryGPU, keepPositive, numAngularVec, numTimeVec)

%% EVALMETRIC: Evaluates the metric and returns the core analysis products
%
%   INPUTS:
%   metric - Metric tensor struct object.
%
%   tryGPU - A flag on whether or not to use GPU computation (0=no, 1=yes)
%
%   keepPositive - A flag on whether or not to return positive values of the energy conditions (0=no, 1=yes).
%
%   numAngularVec - Number of equally spaced spatial vectors to evaluate
%
%   numTimeVec - Number of equally spaced temporal shells to evaluate
%
%
%   OUTPUTS:
%   output - Struct which packages the metric, energy tensors, energy
%   conditions, and scalars.
%

%%

% Handle default input arguments
if nargin < 2
    tryGPU = 0;
end
if nargin < 3
    keepPositive = 1;
end
if nargin < 4
    numAngularVec = 100;
end
if nargin < 5
    numTimeVec = 10;
end

% Metric output
output.metric = metric;

% Energy tensor outputs
output.energyTensor = getEnergyTensor(metric, tryGPU);
output.energyTensorEulerian = doFrameTransfer(metric, output.energyTensor, "Eulerian", tryGPU);

% Energy condition outputs
[output.null] = getEnergyConditions(output.energyTensor, metric, "Null", numAngularVec, numTimeVec, 0, tryGPU);
[output.weak] = getEnergyConditions(output.energyTensor, metric, "Weak", numAngularVec, numTimeVec, 0, tryGPU);
[output.strong] = getEnergyConditions(output.energyTensor, metric, "Strong", numAngularVec, numTimeVec, 0, tryGPU);
[output.dominant] = getEnergyConditions(output.energyTensor, metric, "Dominant", numAngularVec, numTimeVec, 0, tryGPU);

if ~keepPositive
    output.null(output.null > 0) = 0;
    output.weak(output.weak > 0) = 0;
    output.strong(output.strong > 0) = 0;
    output.dominant(output.dominant > 0) = 0;
end

% Scalar outputs
[output.expansion, output.shear, output.vorticity] = getScalars(metric);

end
