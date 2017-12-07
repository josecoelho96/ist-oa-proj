function [ results ] = sensor_validation( A, b, x, threshold, k, s)
results = zeros(1, k);
for i=1:s
    %norm inf must be lower than the threshold for all reliable sensors
    results(i) = norm(b(:,:,i) - A(:,:, i)*x, inf) <= threshold;
end
end