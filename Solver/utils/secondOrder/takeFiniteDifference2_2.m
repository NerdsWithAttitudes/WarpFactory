function [B] = takeFiniteDifference2_2(A,k1,k2,delta)
%TAKEFINITEDIFFERENCE takes the partitial derivatives of the 4-D metric A in the k1 and k2 coordinate direction
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

if s(k1)>=3 && s(k2)>=3
    if k1==k2
        switch k1
            case 1
            B(2:end-1,:,:,:) = (A(3:end,:,:,:)-2*A(2:end-1,:,:,:)+A(1:end-2,:,:,:))/(delta(k1)^2);
            B(1,:,:,:) = B(2,:,:,:);
            B(end,:,:,:) = B(end-1,:,:,:);
        case 2
            B(:,2:end-1,:,:) = (A(:,3:end,:,:)-2*A(:,2:end-1,:,:)+A(:,1:end-2,:,:))/(delta(k1)^2);
            B(:,1,:,:) = B(:,2,:,:);
            B(:,end,:,:) = B(:,end-1,:,:);
        case 3
            B(:,:,2:end-1,:) = (A(:,:,3:end,:)-2*A(:,:,2:end-1,:)+A(:,:,1:end-2,:))/(delta(k1)^2);
            B(:,:,1,:) = B(:,:,2,:);
            B(:,:,end,:) = B(:,:,end-1,:);
        case 4
            B(:,:,:,2:end-1) = (A(:,:,:,3:end)-2*A(:,:,:,2:end-1)+A(:,:,:,1:end-2))/(delta(k1)^2);
            B(:,:,:,1) = B(:,:,:,2);
            B(:,:,:,end) = B(:,:,:,end-1);
        end
    else
        kL = max(k1,k2);
        kS = min(k1,k2);
    
        x1 = 3:s(kS);
        x0 = 2:s(kS)-1;
        x_1 = 1:s(kS)-2;
        
        y1 = 3:s(kL);
        y0 = 2:s(kL)-1;
        y_1 = 1:s(kL)-2;
    
        switch kS
            case 1
                switch kL
                    case 2 %partial t / partial x 
                        B(x0,y0,:,:) = 1/(2^2*delta(kL)*delta(kS)) * ( ...
                            A(x_1,y_1,:,:) - A(x_1,y1,:,:) - A(x1,y_1,:,:) + A(x1,y1,:,:) );
    
                    case 3 %partial t / partial y 
                        B(x0,:,y0,:) = 1/(2^2*delta(kL)*delta(kS)) * ( ...
                            A(x_1,:,y_1,:) - A(x_1,:,y1,:) - A(x1,:,y_1,:) + A(x1,:,y1,:) );

                    case 4 %partial t / partial z
                        B(x0,:,:,y0) = 1/(2^2*delta(kL)*delta(kS)) * ( ...
                            A(x_1,:,:,y_1) - A(x_1,:,:,y1) - A(x1,:,:,y_1) + A(x1,:,:,y1) );
    
                end
            case 2
                switch kL
                    case 3 %partial x / partial y
                        B(:,x0,y0,:) = 1/(2^2*delta(kL)*delta(kS)) * ( ...
                            A(:,x_1,y_1,:) - A(:,x_1,y1,:) - A(:,x1,y_1,:) + A(:,x1,y1,:) );
    
                    case 4 %partial x / partial z
                        B(:,x0,:,y0) = 1/(2^2*delta(kL)*delta(kS)) * ( ...
                            A(:,x_1,:,y_1) - A(:,x_1,:,y1) - A(:,x1,:,y_1) + A(:,x1,:,y1) );
                end
            case 3
                switch kL
                    case 4 %partial y / partial z
                        B(:,:,x0,y0) = 1/(2^2*delta(kL)*delta(kS)) * ( ...
                            A(:,:,x_1,y_1) - A(:,:,x_1,y1) - A(:,:,x1,y_1) + A(:,:,x1,y1) );
                end
        end
    end
end



end

