function [verified] = verifyTensor(inputTensor, suppressMsgs)
%VERIFYMETRIC Verifies the metric tensor and stress energy tensor structs

% Handle input arguments
if nargin < 2
    suppressMsgs = 0;
end

% Get user's backtrace setting
userBacktraceState = warning('query','backtrace').state;
warning('off','backtrace')

verified = 1;

% Check if type field exists
if isfield(inputTensor,'type')
    % Check tensor type
    if strcmpi(inputTensor.type, "Metric")
    
        % Metric tensor
        dispMessage("type: Metric",suppressMsgs);
        
    elseif strcmpi(inputTensor.type, "Stress-Energy")
    
        % Stress-Energy Tensor
        dispMessage("Type: Stress-Energy",suppressMsgs);
    
    elseif ~isfield(inputTensor,'type')
        warning('Tensor type field does not exits. Must be either "Metric" or "Stress-Energy"')
        verified = 0;
    else
        warning("Unknown type")
        verified = 0;
    end


    % Check other properties
    % Tensor
    if isfield(inputTensor,'tensor')
        if isa(inputTensor.tensor,'cell') && size(inputTensor.tensor,1) == 4 && size(inputTensor.tensor,2) == 4 && length(size(inputTensor.tensor{1,1})) == 4
            dispMessage("tensor: Verified",suppressMsgs);
        else
            warning("Tensor is not formatted correctly. Tensor must be a 4x4 cell array of 4D values.")
            verified = 0;
        end
        
    else
        warning("tensor: Empty");
        verified = 0;
    end
    % Coords
    if isfield(inputTensor,'coords')
        if strcmpi(string(inputTensor.coords), "cartesian")
            dispMessage("coords: " + string(inputTensor.coords),suppressMsgs);
        else
            warning("Non-cartesian coordinates are not supported at this time. Set .coords to 'cartesian'.")
        end
    else
        warning("coords: Empty");
        verified = 0;
    end
    % Index
    if isfield(inputTensor,'index')
        if strcmpi(string(inputTensor.index), "contravariant") || ...
           strcmpi(string(inputTensor.index), "covariant") || ...
           strcmpi(string(inputTensor.index), "mixedupdown") || ...
           strcmpi(string(inputTensor.index), "mixeddownup")
            dispMessage("index: " + string(inputTensor.index),suppressMsgs);
        else
            warning("Unknown index");
            verified = 0;
        end
    else
        warning("index: Empty");
        verified = 0;
    end
else
    warning('Tensor type does not exists. Must be Either "Metric" or "Stress-Energy"')
    verified = 0;
end

% Reset user's backtrace setting
warning(userBacktraceState,'backtrace')

% Function for displaying messages
function dispMessage(msg,sM)
    if ~sM
        fprintf(msg+"\n")
    end
end
end

