close all;
clearvars;

k = 16; % Number of sensors
m = 4; % Size of observation vectors b
n = 20; % Size of unknown vector x
reliable_sensors_list = [6 8 10 12 14]; % Number of consistent sensors
restriction_delta = 10^-6;
threshold_value = 10^-4;

s = 14;
    
reliable_sensors = [ones(1, s) zeros(1, k-s)];

%preallocations
bi = zeros(m, 1, k);

% unknown vector is modeled as x0 ~ N(0, n^(-1/2)In)
x0 = mvnrnd(zeros(1, n), n^(-1)*eye(n))';

% Entries of matrix A are drawn independently from N(0, 1)
Ai = randn(m, n, k);

% generate sensor observation data b
% consistent observations
for i=1:s
    bi(:, :, i) = Ai(:, : ,i)*x0;
end

% unreliable sensors
for i=s+1:k
    bi(:, : , i) = mvnrnd(zeros(1, m), eye(m))';
end
        
% CVX problem solving
cvx_begin quiet
    variable x(n, k)
    variable L(n, 1)
    % define cost function
    for i=1:k
        % f(i) = (x(:,i)-L)'*(x(:,i)-L);
        f(i) = norm(x(:,i)-L);
    end
    minimize(sum(f))
    subject to
    for i=1:k
        norm(bi(:,:,i) - Ai(:,:,i)*x(:,i)) <= restriction_delta;
        % (bi(:,:,i)-Ai(:,:,i)*x(:,i))'*(bi(:,:,i)-Ai(:,:,i)*x(:,i)) <= restriction_delta;
    end
cvx_end
        
% check for reliable sensors
method_reliable_sensors = zeros(1, k);
for i=1:s
    %norm inf must be lower than the threshold for all reliable sensors
    method_reliable_sensors(i) = norm(bi(:,:,i) - Ai(:,:, i)*L, inf) <= threshold_value;
end
results = isequal(reliable_sensors, method_reliable_sensors);

disp(reliable_sensors)
disp(method_reliable_sensors)