function [R_munu] = ricciTMem2(gu,gl,delta)
%RICCI calculates the Ricci tensor for a given metric tensor and its
%inverse

% Form from (https://en.wikipedia.org/wiki/Ricci_curvature) introduction section
s = size(gl{1,1});


% Ricci tensor
R_munu = cell(4,4);


%% Precalculate metric derivatives for speed
diff_1_gl = cell(4,4,4);

for i = 1:4
    for j = i:4
        for k = 1:4
            diff_1_gl{i,j,k} = takeFiniteDifference1_2(gl{i,j},k,delta);
        end
    end
end

% Assign symmetric values
for k = 1:4
    diff_1_gl{2,1,k} = diff_1_gl{1,2,k}; 
    diff_1_gl{3,1,k} = diff_1_gl{1,3,k}; 
    diff_1_gl{3,2,k} = diff_1_gl{2,3,k}; 
    diff_1_gl{4,1,k} = diff_1_gl{1,4,k}; 
    diff_1_gl{4,2,k} = diff_1_gl{2,4,k}; 
    diff_1_gl{4,3,k} = diff_1_gl{3,4,k}; 
end

%% Construct Ricci tensor
for i = 1:4
    for j = i:4
        %Begin loop for each i and j of the Ricci Tensor
        if isgpuarray(gu)
            R_munu{i,j} = zeros(s,'gpuArray');
        else
            R_munu{i,j} = zeros(s);
        end

        for a = 1:4
            for b = 1:4
                % First term
                %R_munu{i,j} = R_munu{i,j}-1/2*(diff_2_gl{i,j,a,b}+diff_2_gl{a,b,i,j}-diff_2_gl{i,b,j,a}-diff_2_gl{j,b,i,a}).*gu{a,b};
                R_munu{i,j} = R_munu{i,j}-1/2*(takeFiniteDifference2_2(gl{i,j},a,b,delta)+takeFiniteDifference2_2(gl{a,b},i,j,delta)-takeFiniteDifference2_2(gl{i,b},j,a,delta)-takeFiniteDifference2_2(gl{j,b},i,a,delta)).*gu{a,b};
                for c = 1:4
                    for d = 1:4
                        % Second term
                        R_munu{i,j} = R_munu{i,j}+1/2*(1/2*diff_1_gl{a,c,i}.*diff_1_gl{b,d,j}+diff_1_gl{i,c,a}.*diff_1_gl{j,d,b}-diff_1_gl{i,c,a}.*diff_1_gl{j,b,d}).*gu{a,b}.*gu{c,d};
                        %R_munu{i,j} = R_munu{i,j}+1/2*(1/2*takeFiniteDifference1(gl{a,c},i,delta).*takeFiniteDifference1(gl{b,d},j,delta)+takeFiniteDifference1(gl{i,c},a,delta).*takeFiniteDifference1(gl{j,d},b,delta)-takeFiniteDifference1(gl{i,c},a,delta).*takeFiniteDifference1(gl{j,b},d,delta)).*gu{a,b}.*gu{c,d};
                        % Third term
                        R_munu{i,j} = R_munu{i,j}-1/4*(diff_1_gl{j,c,i}+diff_1_gl{i,c,j}-diff_1_gl{i,j,c}).*(2*diff_1_gl{b,d,a}-diff_1_gl{a,b,d}).*gu{a,b}.*gu{c,d};
                        %R_munu{i,j} = R_munu{i,j}-1/4*(takeFiniteDifference1(gl{j,c},i,delta)+takeFiniteDifference1(gl{i,c},j,delta)-takeFiniteDifference1(gl{i,j},c,delta)).*(2*takeFiniteDifference1(gl{b,d},a,delta)-takeFiniteDifference1(gl{a,b},d,delta)).*gu{a,b}.*gu{c,d};
                    end
                end
            end
        end

    end
end
% Assign symmetric values
R_munu{2,1} = R_munu{1,2}; 
R_munu{3,1} = R_munu{1,3}; 
R_munu{3,2} = R_munu{2,3}; 
R_munu{4,1} = R_munu{1,4}; 
R_munu{4,2} = R_munu{2,4}; 
R_munu{4,3} = R_munu{3,4};

end

