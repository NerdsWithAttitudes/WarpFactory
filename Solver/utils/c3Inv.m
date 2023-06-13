function [invCellArray] = c3Inv(cellArray)
%CINV Finds the inverse of a 3x3 cell array

[h,w] = size(cellArray);
assert(h==3 && w==3, 'Cell array is not 3x3')

r = cellArray;
det = (r{1,1}.*r{2,2}.*r{3,3} - r{1,1}.*r{2,3}.*r{3,2} - r{1,2}.*r{2,1}.*r{3,3} + r{1,2}.*r{2,3}.*r{3,1} + r{1,3}.*r{2,1}.*r{3,2} - r{1,3}.*r{2,2}.*r{3,1});
invCellArray = {1./det.*(r{2,2}.*r{3,3} - r{2,3}.*r{3,2}), 1./det.*(r{1,3}.*r{3,2} - r{1,2}.*r{3,3}), 1./det.*(r{1,2}.*r{2,3} - r{1,3}.*r{2,2}); 1./det.*(r{2,3}.*r{3,1} - r{2,1}.*r{3,3}), 1./det.*(r{1,1}.*r{3,3} - r{1,3}.*r{3,1}), 1./det.*(r{1,3}.*r{2,1} - r{1,1}.*r{2,3}); 1./det.*(r{2,1}.*r{3,2} - r{2,2}.*r{3,1}), 1./det.*(r{1,2}.*r{3,1} - r{1,1}.*r{3,2}), 1./det.*(r{1,1}.*r{2,2} - r{1,2}.*r{2,1})};
end
