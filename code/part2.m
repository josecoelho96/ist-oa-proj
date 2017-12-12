close all;
clearvars;

k = 16; % Number of sensors
m = 4; % Size of observation vectors b
n = 20; % Size of unknown vector x
reliable_sensors_list = [6 8 10 12 14]; % Number of consistent sensors
delta = 10^-6;
threshold = 10^-4;
MCexperiments = 10;

% save outputs
results = zeros(length(reliable_sensors_list), 1);

for s_index=1:length(reliable_sensors_list)
    
    s = reliable_sensors_list(s_index);
    reliable_sensors = [ones(1, s) zeros(1, k-s)];
    
    for j=1:MCexperiments
    
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
            f = 0;
            for i=1:k
                f = f + (x(:,i)-L)'*(x(:,i)-L);
            end
            minimize(f)
            subject to
            for i=1:k
                (bi(:,:,i)-Ai(:,:,i)*x(:,i))'*(bi(:,:,i)-Ai(:,:,i)*x(:,i)) <= delta;
            end
        cvx_end
        
        % check for reliable sensors
        method_reliable_sensors = zeros(1, k);
        for i=1:s
            %norm inf must be lower than the threshold for all reliable sensors
            method_reliable_sensors(i) = norm(bi(:,:,i) - Ai(:,:, i)*L, inf) <= threshold;
        end
        
        results(s_index) = results(s_index) + isequal(reliable_sensors, method_reliable_sensors);
    end
end

% transform results into %
results = (results/MCexperiments) * 100;