function c = trilinearInterp(F,x)
%ref: https://en.wikipedia.org/wiki/Trilinear_interpolation

x = x+10^(-8);
xd = (x(1)-floor(x(1)))/(ceil(x(1))-floor(x(1)));
yd = (x(2)-floor(x(2)))/(ceil(x(2))-floor(x(2))); 
zd = (x(3)-floor(x(3)))/(ceil(x(3))-floor(x(3)));

c00 = F(floor(x(1)),floor(x(2)),floor(x(3)))*(1-xd)+F(ceil(x(1)),floor(x(2)),floor(x(3)))*xd;
c01 = F(floor(x(1)),floor(x(2)),ceil(x(3)))*(1-xd)+F(ceil(x(1)),floor(x(2)),ceil(x(3)))*xd;
c10 = F(floor(x(1)),ceil(x(2)),floor(x(3)))*(1-xd)+F(ceil(x(1)),ceil(x(2)),floor(x(3)))*xd;
c11 = F(floor(x(1)),ceil(x(2)),ceil(x(3)))*(1-xd)+F(ceil(x(1)),ceil(x(2)),ceil(x(3)))*xd;

c0 = c00*(1-yd)+c10*yd;
c1 = c01*(1-yd)+c11*yd;

c = c0*(1-zd)+c1*zd;

end