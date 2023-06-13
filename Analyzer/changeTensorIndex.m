function [outputTensor] = changeTensorIndex(inputTensor, index, metricTensor)

%% CHANGETENSORINDEX: Changes a tensor's index
%
%   INPUTS:
%   inputTensor - Tensor struct to change the index of
%
%   index - Index to change the inputTensor to such as 'covariant',
%   'contravariant', 'mixedupdown', 'mixeddownup'
%
%   metricTensor - Metric struct
%
%
%   OUTPUTS:
%   outputTensor - Tensor struct in the provided index

%%

% Handle default input arguments
if nargin < 3
    if ~strcmpi(inputTensor.type, "metric")
        error("metricTensor is needed as third input when changing index of non-metric tensors.")
    end
else
    if strcmpi(metricTensor.index, "mixedupdown") || strcmpi(metricTensor.index, "mixeddownup")
        error("Metric tensor cannot be used in mixed index.")
    end
end

% Check for if the index transformation exists
if ~(strcmpi(index, "mixedupdown") || strcmpi(index, "mixeddownup") || strcmpi(index, "covariant") || strcmpi(index, "contravariant"))
    error('Transformation selected is not allowed, use either: "covariant", "contravariant", "mixedupdown", "mixeddownup"')
end

%% Transformations
outputTensor = inputTensor;
if strcmpi(inputTensor.type, "metric")
    if (strcmpi(inputTensor.index, "covariant") && strcmpi(index, "contravariant")) || (strcmpi(inputTensor.index, "contravariant") && strcmpi(index, "covariant"))
        outputTensor.tensor = c4Inv(inputTensor.tensor);
    elseif strcmpi(inputTensor.index, "mixedupdown") || strcmpi(inputTensor.index, "mixeddownup")
        error("Input tensor is a Metric tensor of mixed index.")
    elseif strcmpi(index, "mixedupdown") || strcmpi(index, "mixeddownup")
        error("Cannot convert a metric tensor to mixed index.")
    end
else
    % Contravariant/covariant
    if (strcmpi(inputTensor.index, "covariant") && strcmpi(index, "contravariant"))
        if strcmpi(metricTensor.index, "covariant")
            metricTensor.tensor = c4Inv(metricTensor.tensor);
            metricTensor.index = "contravariant";
        end
        outputTensor.tensor = flipIndex(inputTensor, metricTensor);
    elseif (strcmpi(inputTensor.index, "contravariant") && strcmpi(index, "covariant"))
        if strcmpi(metricTensor.index, "contravariant")
            metricTensor.tensor = c4Inv(metricTensor.tensor);
            metricTensor.index = "covariant";
        end
        outputTensor.tensor = flipIndex(inputTensor, metricTensor);
        % To mixed
    elseif strcmpi(inputTensor.index, "contravariant") && strcmpi(index, "mixedupdown")
        if strcmpi(metricTensor.index, "contravariant")
            metricTensor.tensor = c4Inv(metricTensor.tensor);
            metricTensor.index = "covariant";
        end
        outputTensor.tensor = mixIndex2(inputTensor, metricTensor);
    elseif strcmpi(inputTensor.index, "contravariant") && strcmpi(index, "mixeddownup")
        if strcmpi(metricTensor.index, "contravariant")
            metricTensor.tensor = c4Inv(metricTensor.tensor);
            metricTensor.index = "covariant";
        end
        outputTensor.tensor = mixIndex1(inputTensor, metricTensor);
    elseif strcmpi(inputTensor.index, "covariant") && strcmpi(index, "mixedupdown")
        if strcmpi(metricTensor.index, "covariant")
            metricTensor.tensor = c4Inv(metricTensor.tensor);
            metricTensor.index = "contravariant";
        end
        outputTensor.tensor = mixIndex1(inputTensor, metricTensor);
    elseif strcmpi(inputTensor.index, "covariant") && strcmpi(index, "mixeddownup")
        if strcmpi(metricTensor.index, "covariant")
            metricTensor.tensor = c4Inv(metricTensor.tensor);
            metricTensor.index = "contravariant";
        end
        outputTensor.tensor = mixIndex2(inputTensor, metricTensor);
        % From mixed
    elseif strcmpi(inputTensor.index, "mixedupdown") && strcmpi(index, "contravariant")
        if strcmpi(metricTensor.index, "covariant")
            metricTensor.tensor = c4Inv(metricTensor.tensor);
            metricTensor.index = "contravariant";
        end
        outputTensor.tensor = mixIndex2(inputTensor, metricTensor);
    elseif strcmpi(inputTensor.index, "mixedupdown") && strcmpi(index, "covariant")
        if strcmpi(metricTensor.index, "contravariant")
            metricTensor.tensor = c4Inv(metricTensor.tensor);
            metricTensor.index = "covariant";
        end
        outputTensor.tensor = mixIndex1(inputTensor, metricTensor);
    elseif strcmpi(inputTensor.index, "mixeddownup") && strcmpi(index, "covariant")
        if strcmpi(metricTensor.index, "contravariant")
            metricTensor.tensor = c4Inv(metricTensor.tensor);
            metricTensor.index = "covariant";
        end
        outputTensor.tensor = mixIndex2(inputTensor, metricTensor);
    elseif strcmpi(inputTensor.index, "mixeddownup") && strcmpi(index, "contravariant")
        if strcmpi(metricTensor.index, "covariant")
            metricTensor.tensor = c4Inv(metricTensor.tensor);
            metricTensor.index = "contravariant";
        end
        outputTensor.tensor = mixIndex1(inputTensor, metricTensor);
    end
end
outputTensor.index = index;

%% Helper functions

% Flip index
    function tempOutputTensor = flipIndex(inputTensor, metricTensor)
        tempOutputTensor = cell(4, 4);
        for i = 1:4
            for j = 1:4
                tempOutputTensor{i, j} = zeros(size(inputTensor.tensor{i, j}));
                for a = 1:4
                    for b = 1:4
                        tempOutputTensor{i, j} = tempOutputTensor{i, j}+inputTensor.tensor{a, b} .* metricTensor.tensor{a, i} .* metricTensor.tensor{b, j};
                    end
                end
            end
        end
    end

% Mix index 1
    function tempOutputTensor = mixIndex1(inputTensor, metricTensor)
        tempOutputTensor = cell(4, 4);
        for i = 1:4
            for j = 1:4
                tempOutputTensor{i, j} = zeros(size(inputTensor.tensor{i, j}));
                for a = 1:4
                    tempOutputTensor{i, j} = tempOutputTensor{i, j}+inputTensor.tensor{a, j} .* metricTensor.tensor{a, i};
                end
            end
        end
    end

% Mix index 2
    function tempOutputTensor = mixIndex2(inputTensor, metricTensor)
        tempOutputTensor = cell(4, 4);
        for i = 1:4
            for j = 1:4
                tempOutputTensor{i, j} = zeros(size(inputTensor.tensor{i, j}));
                for a = 1:4
                    tempOutputTensor{i, j} = tempOutputTensor{i, j}+inputTensor.tensor{i, a} .* metricTensor.tensor{a, j};
                end
            end
        end
    end

end
