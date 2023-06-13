function [expansionScalar, shearScalar, vorticityScalar] = getScalars(metric)

%% GETSCALARS: Gets the scalars of the input metric
%
%   INPUTS:
%   metric - Metric tensor struct object.
%
%
%   OUTPUTS:
%   expansionScalar - Expansion scalar for all points in spacetime
%   shearScalar - Shear scalar for all points in spacetime
%   vorticityScalar - Vorticity scalar for all points in spacetime
%

%%

% Convert metric tensor to array
arrayMetricTensor(1, 1, 1, 1, :, :) = metric.tensor;
arrayMetricTensor = cell2mat(arrayMetricTensor);

% Get alpha and beta from three plus one decomposer
[alpha, ~, ~, betaUp, ~] = threePlusOneDecomposer(metric);

% Make array of beta at all points in spacetime
arrayBeta(1, 1, 1, 1, :) = betaUp;
arrayBeta = cell2mat(arrayBeta);

% Preallocate u for speed
s = size(metric.tensor{1,1});
uUp = zeros(s(1), s(2), s(3), s(4), 4);
uDown = zeros(s(1), s(2), s(3), s(4), 4);

% Iterate through all points to create
for t = 1:s(1)
    for i = 1:s(2)
        for j = 1:s(3)
            for k = 1:s(4)
                uUp(t, i, j, k, :) = 1 ./ alpha(t, i, j, k) .* [1, -arrayBeta(t, i, j, k, 1), -arrayBeta(t, i, j, k, 2), -arrayBeta(t, i, j, k, 3)];
                uDown(t, i, j, k, :) = squeeze(arrayMetricTensor(t, i, j, k, :, :)) * squeeze(uUp(t, i, j, k, :));
            end
        end
    end
end

% Convert back to cell array
uUpCell = cell(1, 4);
uDownCell = cell(1, 4);
for i = 1:4
    uUpCell{i} = uUp(:, :, :, :, i);
    uDownCell{i} = uDown(:, :, :, :, i);
end

% Calculate covariant derivate
delU = cell(4, 4);
% Make sure metric is covariant
metric = changeTensorIndex(metric, 'covariant');
for i = 1:4
    for j = 1:4
        delU{i, j} = covDiv(metric.tensor, c4Inv(metric.tensor), uUpCell, uDownCell, i, j, [1, 1, 1, 1], 0);
    end
end

% Calcualte projection tensor
P_mix = cell(4, 4);
P = cell(4, 4);
for i = 1:4
    for j = 1:4
        if i == j
            kDelta = 1;
        else
            kDelta = 0;
        end
        % This is correct
        P_mix{i, j} = kDelta+uUpCell{i} .* uDownCell{j};
        P{i, j} = metric.tensor{i, j}+uDownCell{i} .* uDownCell{j};
    end
end

% Define theta tensor
theta.index = "covariant";
theta.type = "tensor";
theta.tensor = cell(4);

% Define omega tensor
omega.index = "covariant";
omega.type = "tensor";
omega.tensor = cell(4);

% Build the tensors
for i = 1:4
    for j = 1:4
        theta.tensor{i, j} = zeros(size(metric.tensor{1, 1}));
        omega.tensor{i, j} = zeros(size(metric.tensor{1, 1}));
        for a = 1:4
            for b = 1:4
                theta.tensor{i, j} = theta.tensor{i, j}+P_mix{a, i} .* P_mix{b, j} .* 1 / 2 .* (delU{a, b}+delU{b, a}); % \nabla_{(\alpha} U_{\beta)}=\frac{1}{2}\left(\nabla_\alpha U_\beta+\nabla_\beta U_\alpha\right)
                omega.tensor{i, j} = omega.tensor{i, j}+P_mix{a, i} .* P_mix{b, j} .* 1 / 2 .* (delU{a, b}-delU{b, a}); % \nabla_{[\alpha} U_{\beta]}=\frac{1}{2}\left(\nabla_\alpha U_\beta-\nabla_\beta U_\alpha\right)
            end
        end
    end
end

% Get the trace of theta to calculate its scalar
thetaTrace = getTrace(theta, metric);

% Calculate omega scalar
omega_up = changeTensorIndex(omega, "contravariant", metric);
omegaTrace = zeros(size(metric.tensor{1, 1}));
for mu = 1:4
    for nu = 1:4
        omegaTrace = omegaTrace+1 / 2 * omega_up.tensor{mu, nu} .* omega.tensor{mu, nu};
    end
end

%% Shear
% Define shear tensor
shear.index = "covariant";
shear.type = "tensor";
shear.tensor = cell(4);
for i = 1:4
    for j = 1:4
        shear.tensor{i, j} = theta.tensor{i, j}-thetaTrace .* 1 / 3 .* P{i, j};
    end
end

% Calculate scalar shear
shear_up = changeTensorIndex(shear, "contravariant", metric);
sigma2 = zeros(size(metric.tensor{1, 1}));
for i = 1:4
    for j = 1:4
        sigma2 = sigma2+1 / 2 * shear.tensor{i, j} .* shear_up.tensor{i, j};
    end
end

% Output results
shearScalar = sigma2;
expansionScalar = thetaTrace;
vorticityScalar = omegaTrace;

end
