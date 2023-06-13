function indexData = getSliceData(plane, sliceCenter, tensor)
    s = size(tensor.tensor{1,1});
    indexData = {1:s(1); 1:s(2); 1:s(3); 1:s(4)};

    indexData{plane(1)} = sliceCenter(1);
    indexData{plane(2)} = sliceCenter(2);
end