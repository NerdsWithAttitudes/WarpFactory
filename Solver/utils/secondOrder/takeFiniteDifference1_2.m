function [B] = takeFiniteDifference1_2(A,k,delta)
%TAKEFINITEDIFFERENCE takes the partitial derivative of the 4-D metric A in the k coordinate direction

if isgpuarray(A)
    delta = gpuArray(delta);

    s = gpuArray(size(A));
    B = zeros(s,'gpuArray');
else
    s = size(A);
    B = zeros(s);
end

if s(k)>=3
    switch k
        case 1
            B(2:end-1,:,:,:) = (A(3:end,:,:,:)-A(1:end-2,:,:,:))/(2*delta(k));
            B(1,:,:,:) = B(2,:,:,:);
            B(end,:,:,:) = B(end-1,:,:,:);
        case 2
            B(:,2:end-1,:,:) = (A(:,3:end,:,:)-A(:,1:end-2,:,:))/(2*delta(k));
            B(:,1,:,:) = B(:,2,:,:);
            B(:,end,:,:) = B(:,end-1,:,:);
        case 3
            B(:,:,2:end-1,:) = (A(:,:,3:end,:)-A(:,:,1:end-2,:))/(2*delta(k));
            B(:,:,1,:) = B(:,:,2,:);
            B(:,:,end,:) = B(:,:,end-1,:);
        case 4
            B(:,:,:,2:end-1) = (A(:,:,:,3:end)-A(:,:,:,1:end-2))/(2*delta(k));
            B(:,:,:,1) = B(:,:,:,2);
            B(:,:,:,end) = B(:,:,:,end-1);
    end
end


end

