function [Vector] = getEvenPointsOnSphere(R,numberOfPoints,useGPU)
%GETEVENPOINTSONSPHERE Summary of this function goes here
%   Detailed explanation goes here

if useGPU
    goldenRatio = gpuArray((1+5^0.5)/2);
    numberOfPoints = gpuArray(numberOfPoints);
    R = gpuArray(R);
    Vector = zeros(3,numberOfPoints,'gpuArray');
else
    goldenRatio = (1+5^0.5)/2;
    Vector = zeros(3,numberOfPoints);
end

for i = 0:numberOfPoints-1
    theta = 2*pi*i/goldenRatio;
    phi = acos(1-2*(i+0.5)/numberOfPoints);

    Vector(1,i+1) = R*cos(theta)*sin(phi);
    Vector(2,i+1) = R*sin(theta)*sin(phi);
    Vector(3,i+1) = R*cos(phi);
end

Vector = real(Vector);
end

