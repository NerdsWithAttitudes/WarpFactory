function arrayTensor = tensorCell2Array(Tensor,tryGPU)
% Converts a 4x4 cells of 4D spactime array into a single array with indexing of (mu,nu,t,x1,x2,x3)
if nargin < 2
    tryGPU = 0;
end

    arrayTensor(1,1,1,1,:,:) = Tensor.tensor;
    if tryGPU
        arrayTensor = cell2matGPU(arrayTensor);
    else
        arrayTensor = cell2mat(arrayTensor);
    end
end