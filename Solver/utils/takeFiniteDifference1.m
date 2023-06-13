function [B] = takeFiniteDifference1(A,k,delta,phiphiFlag)
%TAKEFINITEDIFFERENCE takes the partial derivative of the 4-D metric A in the k coordinate direction

if isgpuarray(A)
    delta = gpuArray(delta);

    s = gpuArray(size(A));
    B = zeros(s,'gpuArray');
else
    s = size(A);
    B = zeros(s);
end

if s(k)>=5
    switch k
        case 1
            B(3:end-2,:,:,:) = (-(A(5:end,:,:,:)-A(1:end-4,:,:,:))+8*(A(4:end-1,:,:,:)-A(2:end-3,:,:,:)))/(12*delta(k));
            B(1,:,:,:) = B(3,:,:,:);
            B(2,:,:,:) = B(3,:,:,:);
            B(end-1,:,:,:) = B(end-2,:,:,:);
            B(end,:,:,:) = B(end-2,:,:,:);
        case 2
            B(:,3:end-2,:,:) = (-(A(:,5:end,:,:)-A(:,1:end-4,:,:))+8*(A(:,4:end-1,:,:)-A(:,2:end-3,:,:)))/(12*delta(k));
            B(:,1,:,:) = B(:,3,:,:);
            B(:,2,:,:) = B(:,3,:,:);
            B(:,end-1,:,:) = B(:,end-2,:,:);
            B(:,end,:,:) = B(:,end-2,:,:);
        case 3
            B(:,:,3:end-2,:) = (-(A(:,:,5:end,:)-A(:,:,1:end-4,:))+8*(A(:,:,4:end-1,:)-A(:,:,2:end-3,:)))/(12*delta(k));
            if phiphiFlag
                B(:,:,1,:) = 2*4;
                B(:,:,2,:) = 2*3;
                B(:,:,end-1,:) = 2*(s(3)-5-1);
                B(:,:,end,:) = 2*(s(3)-5);
            else
                B(:,:,end-1,:) = B(:,:,end-2,:);
                B(:,:,end,:) = B(:,:,end-2,:);
                B(:,:,1,:) = B(:,:,3,:);
                B(:,:,2,:) = B(:,:,3,:);
            end
        case 4
            B(:,:,:,3:end-2) = (-(A(:,:,:,5:end)-A(:,:,:,1:end-4))+8*(A(:,:,:,4:end-1)-A(:,:,:,2:end-3)))/(12*delta(k));
            B(:,:,:,1) = B(:,:,:,3);
            B(:,:,:,2) = B(:,:,:,3);
            B(:,:,:,end-1) = B(:,:,:,end-2);
            B(:,:,:,end) = B(:,:,:,end-2);
    end
end


end

