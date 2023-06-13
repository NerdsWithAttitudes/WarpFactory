function cdVec = covDiv(gl, gu, vecU, vecD, idxDiv, idxVec, delta, stairSel)

% Precalculate metric derivatives
diff_1_gl = cell(4, 4, 4);
s = size(gl{1, 1});
for i = 1:4
    for j = 1:4
        if i == 2 && j == 2 && s(2) == 1
            phiphiFlag = 1;
        else
            phiphiFlag = 0;
        end
        for k = 1:4
            diff_1_gl{i, j, k} = takeFiniteDifference1(gl{i, j}, k, delta, phiphiFlag);
        end
    end
end

% Covariant derivative of covariant vector
if stairSel == 0 
    % Build gradient operated vector
    cdVec = takeFiniteDifference1(vecD{idxVec}, idxDiv, delta, 0);
    for i = 1:4
        Gamma = getChristoffelSym(gu, diff_1_gl, i, idxVec, idxDiv);
        cdVec = cdVec-Gamma .* vecD{i};
    end
elseif stairSel == 1 %covariant dereivative of contravariant vector
    % Build gradient operated vector
    cdVec = takeFiniteDifference1(vecU{idxVec}, idxDiv, delta, 0);
    for i = 1:4
        Gamma = getChristoffelSym(gu, diff_1_gl, idxVec, idxDiv, i);
        cdVec = cdVec+Gamma .* vecU{i};
    end
else
    error('Invalid variance selected')
end


end