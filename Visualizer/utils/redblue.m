function [returnedMap] = redblue(value, gradientNum)

if nargin < 2
    gradientNum = 1024;
end

minValue = min(value,[],'all');
maxValue = max(value,[],'all');

if ~(minValue <=0 && maxValue >=0)
    if minValue > 0 && maxValue > 0
        returnedMap = cat(1,linspace(1,0,round(gradientNum)),linspace(1,0,round(gradientNum)),ones(1,round(gradientNum)))';
        return
    end
    if minValue < 0 && maxValue < 0
        returnedMap = cat(1,ones(1,round(gradientNum)),linspace(0,1,round(gradientNum)),linspace(0,1,round(gradientNum)))';
        return
    end
end

if minValue == 0 && maxValue == 0
    returnedMap = [1 1 1];
    return
end

centerVal = abs(maxValue)/(abs(maxValue)+abs(minValue));
returnedMap = cat(1,ones(1,round((1-centerVal)*gradientNum)),linspace(0,1,round((1-centerVal)*gradientNum)),linspace(0,1,round((1-centerVal)*gradientNum)))';
returnedMap = [returnedMap; cat(1,linspace(1,0,round(centerVal*gradientNum)),linspace(1,0,round(centerVal*gradientNum)),ones(1,round(centerVal*gradientNum)))'];

end

