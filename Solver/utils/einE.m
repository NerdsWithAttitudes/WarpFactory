function [enDen] = einE(E,gu)
%EIND Calculates the Energy Tensor from the Einstien Tensor
% Returns value in Joules/m^3

enDen_ = cell(4,4);
for mu = 1:4
    for nu = 1:4
        enDen_{mu,nu} = c^4./(8.*pi.*G).*E{mu,nu};
    end
end

% Turn into contravarient form
enDen = cell(4,4);
for mu = 1:4
    for nu = 1:4
        enDen{mu,nu} = 0;
        for alpha = 1:4
            for beta = 1:4
                enDen{mu,nu} = enDen{mu,nu} + enDen_{alpha,beta}.*gu{alpha,mu}.*gu{beta,nu};
            end
        end
    end
end

end

