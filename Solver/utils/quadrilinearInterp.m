function c = quadrilinearInterp(F,x)
% runs a weighted average between two timeslices after a spatial trilinear
% interpolation in both slices.
t1 = floor(x(1));
t2 = ceil(x(1));
t = x(1);

if t1 == t2
    c = trilinearInterp(squeeze(F(t1,:,:,:)),x(2:end));
else
    c1 = trilinearInterp(squeeze(F(t1,:,:,:)),x(2:end));
    c2 = trilinearInterp(squeeze(F(t2,:,:,:)),x(2:end));
    
    c = (c1*(t2-t)+c2*(t-t1))/(t2-t1);
end

end