function [xlabelName, ylabelName] = labelCartesianAxis(plane)
    labels = ["t", "x", "y", "z"];
    shownPlanes = sort(setdiff(1:4, plane));

    xlabelName = labels(shownPlanes(1));
    ylabelName = labels(shownPlanes(2));
end
