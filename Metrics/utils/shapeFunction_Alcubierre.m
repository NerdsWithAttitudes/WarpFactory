function [f] = shapeFunction_Alcubierre(r,R,sigma)
    f = (tanh(sigma*(R + r)) + tanh(sigma*(R - r)))/(2*tanh(R*sigma));
end

