function VecField = generateUniformField(type,numAngularVec,numTimeVec,tryGPU)

if ~(strcmpi(type,"nulllike") || strcmpi(type,"timelike"))
    error('Vector field type not generated, use either: "nulllike", "timelike"')
end

if strcmpi(type,"timelike")
    % generate timelike vectors c^2t^2 > r^2, 
    bb = linspace(0,1,numTimeVec);
    VecField = ones(4,numAngularVec,numTimeVec);
    for jj = 1:numTimeVec
        % build vector field in cartesian coordinates
        VecField(:,:,jj) = [ones(1,numAngularVec); getEvenPointsOnSphere(1-bb(jj),numAngularVec,1)];
        VecField(:,:,jj) = VecField(:,:,jj)./(VecField(1,:,jj).^2+VecField(2,:,jj).^2+VecField(3,:,jj).^2+VecField(4,:,jj).^2).^0.5;
    end
 

elseif strcmpi(type,"nulllike")
    % build vector field in cartesian coordinates
    getEvenPointsOnSphere(1,numAngularVec,1);
    VecField = ones(4,numAngularVec);
    VecField(2:end,:) = getEvenPointsOnSphere(1,numAngularVec,1);
    VecField = VecField./(VecField(1,:).^2+VecField(2,:).^2+VecField(3,:).^2+VecField(4,:).^2).^0.5;
end

% Declare variables to be determined in eval of E conditions
if tryGPU
    VecField = gpuArray(VecField);
end

end