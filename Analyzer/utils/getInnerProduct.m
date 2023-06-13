function innerprod = getInnerProduct(vecA,vecB,Metric)
% takes the innerproduct of two vector fields with their metric (either
% array form or struct form

if ~verifyTensor(Metric,1)
    error("Metric is not verified. Please verify metric using verifyTensor(metric).")
end

s = size(Metric.tensor{1,1});
innerprod = zeros(s);

if ~strcmpi(vecA.index,vecB.index)
    for mu = 1:4
        for nu = 1:4
            innerprod = innerprod + vecA.field{mu}.*vecB.field{nu};
        end
    end

elseif strcmpi(vecA.index,vecB.index)
    if strcmpi(vecA.index, Metric.index)
        Metric.tensor = c4Inv(Metric.tensor); %flip index 
    end
    for mu = 1:4
        for nu = 1:4
            innerprod = innerprod + vecA.field{mu}.*vecB.field{nu}.*Metric.tensor{mu,nu};
        end
    end
end

      
end