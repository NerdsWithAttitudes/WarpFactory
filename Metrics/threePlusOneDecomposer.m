function [alpha, betaDown, gammaDown, betaUp, gammaUp] = threePlusOneDecomposer(metric)
%% THREEPLUSONEDECOMPOSER: Finds 3+1 terms from the metric tensor 
%
%   INPUTS:
%   metric - metric struct object. 
%
%   OUPUTS: 
%   alpha - 4D array. Lapse rate.
%
%   betaDown - 1x3 cell of 4D arrays. Shift vectors.
%
%   gamma - 3x3 cell of 4D arrays. Spatial terms.

%%

% Check that the metric is covariant and change index if not
metric = changeTensorIndex(metric, "covariant");

% Covariant shift vector maps to the covariant tensor terms g_0i
betaDown = {metric.tensor{1,2}; metric.tensor{1,3}; metric.tensor{1,4}};

% Covariant gamma maps to the covariant tensor terms g_ij
gammaDown = {metric.tensor{2,2} metric.tensor{2,3}, metric.tensor{2,4}; ...
        metric.tensor{3,2}, metric.tensor{3,3}, metric.tensor{3,4}; ...
        metric.tensor{4,2}, metric.tensor{4,3}, metric.tensor{4,4}};

% Transform gamma to contravariant
gammaUp = c3Inv(gammaDown);

% Find the world gridSize
s = size(metric.tensor{1,1});

% Transform beta to contravariant
betaUp = cell(1,3);
for i = 1:3
    betaUp{i} = zeros(s); 
    for j = 1:3
        betaUp{i} = betaUp{i} + gammaUp{i,j}.*betaDown{j};
    end
end

% Find lapse using beta and g_00
alpha = sqrt(betaUp{1}.*betaDown{1} + betaUp{2}.*betaDown{2} + betaUp{3}.*betaDown{3} - metric.tensor{1,1});

end

