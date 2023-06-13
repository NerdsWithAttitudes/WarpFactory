function [Trace] = getTrace(tensor,metric)
%GETTRACE take the trace of a tensor

% Verify metric
if ~verifyTensor(metric,1)
    error("Metric is not verified. Please verify metric using verifyTensor(metric).")
end

Trace = zeros(size(metric.tensor{1,1}));

if strcmpi(tensor.index,metric.index)
    metric.tensor = c4Inv(metric.tensor);
end

for a = 1:4
    for b = 1:4
        Trace = Trace + metric.tensor{a,b}.*tensor.tensor{a,b};
    end
end
           

end

