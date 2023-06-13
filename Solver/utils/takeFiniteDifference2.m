function [B] = takeFiniteDifference2(A,k1,k2,delta,phiphiFlag)
%TAKEFINITEDIFFERENCE takes the partial derivatives of the 4-D metric A in the k1 and k2 coordinate direction
% Assumes that boundary terms are zero from constant first order
% derivatives


if isgpuarray(A)
    delta = gpuArray(delta);

    s = gpuArray(size(A));
    B = zeros(s,'gpuArray');
    
else
    s = size(A);
    B = zeros(s);
end

if s(k1)>=5 && s(k2)>=5
    if k1==k2
        switch k1
            case 1
                B(3:end-2,:,:,:) = (-(A(5:end,:,:,:)+A(1:end-4,:,:,:))+16*(A(4:end-1,:,:,:)+A(2:end-3,:,:,:))-30*A(3:end-2,:,:,:))/(12*delta(k1)^2);
                B(1,:,:,:) = B(3,:,:,:);
                B(2,:,:,:) = B(3,:,:,:);
                B(end-1,:,:,:) = B(end-2,:,:,:);
                B(end,:,:,:) = B(end-2,:,:,:);
            case 2
                B(:,3:end-2,:,:) = (-(A(:,5:end,:,:)+A(:,1:end-4,:,:))+16*(A(:,4:end-1,:,:)+A(:,2:end-3,:,:))-30*A(:,3:end-2,:,:))/(12*delta(k1)^2);
                B(:,1,:,:) = B(:,3,:,:);
                B(:,2,:,:) = B(:,3,:,:);
                B(:,end-1,:,:) = B(:,end-2,:,:);
                B(:,end,:,:) = B(:,end-2,:,:);
            case 3
                B(:,:,3:end-2,:) = (-(A(:,:,5:end,:)+A(:,:,1:end-4,:))+16*(A(:,:,4:end-1,:)+A(:,:,2:end-3,:))-30*A(:,:,3:end-2,:))/(12*delta(k1)^2);
                if phiphiFlag
                    B(:,:,1,:) = -2;
                    B(:,:,2,:) = -2;
                    B(:,:,end-1,:) = 2;
                    B(:,:,end,:) = 2;
                else
                    B(:,:,end-1,:) = B(:,:,end-2,:);
                    B(:,:,end,:) = B(:,:,end-2,:);
                    B(:,:,1,:) = B(:,:,3,:);
                    B(:,:,2,:) = B(:,:,3,:);
                end
            case 4
                B(:,:,:,3:end-2) = (-(A(:,:,:,5:end)+A(:,:,:,1:end-4))+16*(A(:,:,:,4:end-1)+A(:,:,:,2:end-3))-30*A(:,:,:,3:end-2))/(12*delta(k1)^2);
                B(:,:,:,1) = B(:,:,:,3);
                B(:,:,:,2) = B(:,:,:,3);
                B(:,:,:,end-1) = B(:,:,:,end-2);
                B(:,:,:,end) = B(:,:,:,end-2);
        end
    else
        kL = max(k1,k2);
        kS = min(k1,k2);
    
        x2 = 5:s(kS);
        x1 = 4:s(kS)-1;
        x0 = 3:s(kS)-2;
        x_1 = 2:s(kS)-3;
        x_2 = 1:s(kS)-4;
        
        y2 = 5:s(kL);
        y1 = 4:s(kL)-1;
        y0 = 3:s(kL)-2;
        y_1 = 2:s(kL)-3;
        y_2 = 1:s(kL)-4;
    
        switch kS
            case 1
                switch kL
                    case 2 %partial 0 / partial 1
                        B(x0,y0,:,:) = 1/(12^2*delta(kL)*delta(kS)) * ( ...
                            -(-(A(x2,y2,:,:)  -A(x_2,y2,:,:) ) +8*(A(x1,y2,:,:)  -A(x_1,y2,:,:) ))...
                            +(-(A(x2,y_2,:,:) -A(x_2,y_2,:,:)) +8*(A(x1,y_2,:,:) -A(x_1,y_2,:,:)))...
                          +8*(-(A(x2,y1,:,:)  -A(x_2,y1,:,:) ) +8*(A(x1,y1,:,:)  -A(x_1,y1,:,:) ))...
                          -8*(-(A(x2,y_1,:,:) -A(x_2,y_1,:,:)) +8*(A(x1,y_1,:,:) -A(x_1,y_1,:,:)))...
                          );
    
                    case 3 %partial 0 / partial 2
                        B(x0,:,y0,:) = 1/(12^2*delta(kL)*delta(kS)) * ( ...
                            -(-(A(x2,:,y2,:)  -A(x_2,:,y2,:) ) +8*(A(x1,:,y2,:)  -A(x_1,:,y2,:) ))...
                            +(-(A(x2,:,y_2,:) -A(x_2,:,y_2,:)) +8*(A(x1,:,y_2,:) -A(x_1,:,y_2,:)))...
                          +8*(-(A(x2,:,y1,:)  -A(x_2,:,y1,:) ) +8*(A(x1,:,y1,:)  -A(x_1,:,y1,:) ))...
                          -8*(-(A(x2,:,y_1,:) -A(x_2,:,y_1,:)) +8*(A(x1,:,y_1,:) -A(x_1,:,y_1,:)))...
                          );
                    case 4 %partial 0 / partial 3
                        B(x0,:,:,y0) = 1/(12^2*delta(kL)*delta(kS)) * ( ...
                            -(-(A(x2,:,:,y2)  -A(x_2,:,:,y2) ) +8*(A(x1,:,:,y2)  -A(x_1,:,:,y2) ))...
                            +(-(A(x2,:,:,y_2) -A(x_2,:,:,y_2)) +8*(A(x1,:,:,y_2) -A(x_1,:,:,y_2)))...
                          +8*(-(A(x2,:,:,y1)  -A(x_2,:,:,y1) ) +8*(A(x1,:,:,y1)  -A(x_1,:,:,y1) ))...
                          -8*(-(A(x2,:,:,y_1) -A(x_2,:,:,y_1)) +8*(A(x1,:,:,y_1) -A(x_1,:,:,y_1)))...
                          );
    
                end
            case 2
                switch kL
                    case 3 %partial 1 / partial 2
                        B(:,x0,y0,:) = 1/(12^2*delta(kL)*delta(kS)) * ( ...
                            -(-(A(:,x2,y2,:)  -A(:,x_2,y2,:) ) +8*(A(:,x1,y2,:)  -A(:,x_1,y2,:) ))...
                            +(-(A(:,x2,y_2,:) -A(:,x_2,y_2,:)) +8*(A(:,x1,y_2,:) -A(:,x_1,y_2,:)))...
                          +8*(-(A(:,x2,y1,:)  -A(:,x_2,y1,:) ) +8*(A(:,x1,y1,:)  -A(:,x_1,y1,:) ))...
                          -8*(-(A(:,x2,y_1,:) -A(:,x_2,y_1,:)) +8*(A(:,x1,y_1,:) -A(:,x_1,y_1,:)))...
                          );
    
                    case 4 %partial 1 / partial 3
                         B(:,x0,:,y0) = 1/(12^2*delta(kL)*delta(kS)) * ( ...
                            -(-(A(:,x2,:,y2)  -A(:,x_2,:,y2) ) +8*(A(:,x1,:,y2)  -A(:,x_1,:,y2) ))...
                            +(-(A(:,x2,:,y_2) -A(:,x_2,:,y_2)) +8*(A(:,x1,:,y_2) -A(:,x_1,:,y_2)))...
                          +8*(-(A(:,x2,:,y1)  -A(:,x_2,:,y1) ) +8*(A(:,x1,:,y1)  -A(:,x_1,:,y1) ))...
                          -8*(-(A(:,x2,:,y_1) -A(:,x_2,:,y_1)) +8*(A(:,x1,:,y_1) -A(:,x_1,:,y_1)))...
                          );
                end
            case 3
                switch kL
                    case 4 %partial 2 / partial 3
                        B(:,:,x0,y0) = 1/(12^2*delta(kL)*delta(kS)) * ( ...
                            -(-(A(:,:,x2,y2)  -A(:,:,x_2,y2) ) +8*(A(:,:,x1,y2)  -A(:,:,x_1,y2) ))...
                            +(-(A(:,:,x2,y_2) -A(:,:,x_2,y_2)) +8*(A(:,:,x1,y_2) -A(:,:,x_1,y_2)))...
                          +8*(-(A(:,:,x2,y1)  -A(:,:,x_2,y1) ) +8*(A(:,:,x1,y1)  -A(:,:,x_1,y1) ))...
                          -8*(-(A(:,:,x2,y_1) -A(:,:,x_2,y_1)) +8*(A(:,:,x1,y_1) -A(:,:,x_1,y_1)))...
                          );
                end
        end
    end
end



end

