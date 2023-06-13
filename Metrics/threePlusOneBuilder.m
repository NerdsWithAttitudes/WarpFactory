function [metricTensor] = threePlusOneBuilder(alpha, beta, gamma)

%% THREEPLUSONEBUILDER: Builds the metric given input 3+1 components of alpha, beta, and gamma
%
%   INPUTS:
%   alpha - (TxXxYxZ) lapse rate map across spacetime
%   beta - {3}x(TxXxYxZ) (covariant assumed) shift vector map across spacetime
%   gamma - {3x3}x(TxXxYxZ) (covariant assumed) spatial term map map across spacetime
%
%
%   OUTPUTS:
%   metricTensor - metric struct
%

%%

% Set spatial components
gamma_up = c3Inv(gamma);

% Find gridSize
s = size(alpha);

% Caluculate beta_i
beta_up = cell(1,3);

for i = 1:3
    beta_up{i} = zeros(s);
    for j = 1:3
        beta_up{i} = beta_up{i} + gamma_up{i, j} .* beta{j};
    end
end

% Create time-time component
metricTensor{1, 1} = -alpha.^2;
for i = 1:3
    metricTensor{1, 1} = metricTensor{1, 1}+beta_up{i} .* beta{i};
end

% Create time-space components
for i = 2:4
    metricTensor{1, i} = beta{i-1};
    metricTensor{i, 1} = metricTensor{1, i};
end

% Create space-space components
for i = 2:4
    for j = 2:4
        metricTensor{i, j} = gamma{i-1, j-1};
    end
end


end
