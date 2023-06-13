function outputValue = legendreRadialInterp(inputArray,r)

% 3rd Order Legendre Polynomial Interpolation
    rScale = 1;

    x0 = floor(r/rScale-1);
    x1 = floor(r/rScale);
    x2 = ceil(r/rScale);
    x3 = ceil(r/rScale+1);

    y0 = inputArray(max(x0,1));
    y1 = inputArray(max(x1,1));
    y2 = inputArray(max(x2,1));
    y3 = inputArray(max(x3,1));
    
    x = r;
    
    x0 = x0*rScale;
    x1 = x1*rScale;
    x2 = x2*rScale;
    x3 = x3*rScale;

    outputValue = y0*(x-x1)*(x-x2)*(x-x3)/((x0-x1)*(x0-x2)*(x0-x3))...
        + y1*(x-x0)*(x-x2)*(x-x3)/((x1-x0)*(x1-x2)*(x1-x3))...
        + y2*(x-x0)*(x-x1)*(x-x3)/((x2-x0)*(x2-x1)*(x2-x3))...
        + y3*(x-x0)*(x-x1)*(x-x2)/((x3-x0)*(x3-x1)*(x3-x2));

end