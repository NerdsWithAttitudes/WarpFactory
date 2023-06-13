function [cellDet] = cDet(cellArray)
%CINV Finds the determinant of a cell array
[h,w] = size(cellArray);
if h==2 && w==2
   cellDet = cellArray{1,1}.*cellArray{2,2}-cellArray{1,2}.*cellArray{2,1};
   return
end
cellDet = 0;
for i = 1:h
    subArray = cellArray;
    subArray(1,:) = [];
    subArray(:,i) = [];
    subDet = cDet(subArray);
    cellDet = cellDet + (2*mod(i,2)-1).*cellArray{1,i}.*subDet;
end
end