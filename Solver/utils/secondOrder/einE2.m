function [enDen] = einE2(E,gu,units)
%EIND Calculates the Energy Tensor from the Einstien Tensor
% Returns value in Joules/m^3
% G = 6.674*10^-11 m^3 kg^-1 s^-2
% c = 2.99792*10^8 m s^-1
G = 6.674*10^-11 * units(3)^2 * units(2) / units(1)^3;
c = 2.99792*10^8 * units(3)   / units(1);

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

