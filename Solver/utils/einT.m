function [E] = einT(R_munu,R,gl)
%EINT Calculates the Eienstien Tensor from the Ricci Tensor, Ricci Scalar,
%and inverse metric tensor

E = cell(4,4);
for mu = 1:4
    for nu = 1:4
        E{mu,nu} = R_munu{mu,nu}-0.5.*gl{mu,nu}.*R;
    end
end

end

