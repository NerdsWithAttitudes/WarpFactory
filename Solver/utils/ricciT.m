function [R_munu] = ricciT(gu,gl,delta)
%RICCI calculates the Ricci tensor for a given metric tensor and its
%inverse

% Form from (https://en.wikipedia.org/wiki/Ricci_curvature) introduction section
s = size(gl{1,1});


% Ricci tensor
R_munu = cell(4,4);


%% Precalculate metric derivatives for speed
diff_1_gl = cell(4,4,4);
diff_2_gl = cell(4,4,4,4);

for i = 1:4
    for j = i:4
        phiphiFlag = 0;
        for k = 1:4
            diff_1_gl{i,j,k} = takeFiniteDifference1(gl{i,j},k,delta,phiphiFlag);
            if k == 1
                diff_1_gl{i,j,k} = 1/c*diff_1_gl{i,j,k};
            end
            for n = k:4
                diff_2_gl{i,j,k,n} = takeFiniteDifference2(gl{i,j},k,n,delta,phiphiFlag);

                if (n == 1 && k ~= 1) || (n ~= 1 && k == 1)
                    diff_2_gl{i,j,k,n} = 1/c*diff_2_gl{i,j,k,n};
                elseif n == 1 && k == 1
                    diff_2_gl{i,j,k,n} = 1/c^2*diff_2_gl{i,j,k,n};
                end

                if k~=n
                    diff_2_gl{i,j,n,k} = diff_2_gl{i,j,k,n};
                end
            end
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
    for n= 1:4
        diff_2_gl{2,1,k,n} = diff_2_gl{1,2,k,n}; 
        diff_2_gl{3,1,k,n} = diff_2_gl{1,3,k,n}; 
        diff_2_gl{3,2,k,n} = diff_2_gl{2,3,k,n}; 
        diff_2_gl{4,1,k,n} = diff_2_gl{1,4,k,n}; 
        diff_2_gl{4,2,k,n} = diff_2_gl{2,4,k,n}; 
        diff_2_gl{4,3,k,n} = diff_2_gl{3,4,k,n}; 
    end
end



%% Construct Ricci tensor
for i = 1:4
    for j = i:4
        % Begin loop for each i and j of the Ricci Tensor
        if isgpuarray(gu{1,1})
            R_munu_temp = zeros(s,'gpuArray');
        else
            R_munu_temp = zeros(s);
        end

        for a = 1:4

            for b = 1:4
                if isgpuarray(gu{1,1})
                    R_munu_temp_2 = zeros(s,'gpuArray');
                else
                    R_munu_temp_2 = zeros(s);
                end

                % First term
                R_munu_temp_2 = R_munu_temp_2 - (diff_2_gl{i,j,a,b}+diff_2_gl{a,b,i,j}-diff_2_gl{i,b,j,a}-diff_2_gl{j,b,i,a});
                
                for r = 1:4
                    if isgpuarray(gu{1,1})
                        R_munu_temp_3 = zeros(s,'gpuArray');
                        R_munu_temp_4 = zeros(s,'gpuArray');
                        R_munu_temp_5 = zeros(s,'gpuArray');
                    else
                        R_munu_temp_3 = zeros(s);
                        R_munu_temp_4 = zeros(s);
                        R_munu_temp_5 = zeros(s);
                    end

                    for d = 1:4
                        % Second term
                        R_munu_temp_3 = R_munu_temp_3 + diff_1_gl{b,d,j}.*gu{r,d};
                        R_munu_temp_4 = R_munu_temp_4 + (diff_1_gl{j,d,b}-diff_1_gl{j,b,d}).*gu{r,d};
                        % Third term
                        R_munu_temp_5 = R_munu_temp_5 - (diff_1_gl{b,d,a}+diff_1_gl{b,d,a}-diff_1_gl{a,b,d}).*gu{r,d};
                    end
                    R_munu_temp_2 = R_munu_temp_2 + R_munu_temp_4.*diff_1_gl{i,r,a} + 1/2*(R_munu_temp_3.*diff_1_gl{a,r,i} + R_munu_temp_5.*(diff_1_gl{j,r,i}+diff_1_gl{i,r,j}-diff_1_gl{j,i,r}));
                end
                R_munu_temp = R_munu_temp + gu{a,b}.*R_munu_temp_2;
            end
        end
        R_munu{i,j} = 1/2*R_munu_temp;
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

