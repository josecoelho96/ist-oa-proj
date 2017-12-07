% Optimization and Algorithms
% Project 2
% 2017/2018 Fall Semester
% Group 23
% José Coelho, 81013
% Miguel de Moura, 78887
% Gonçalo Pereira, 81602

% clear workspace and close all figures
close all;
clear;

% Setting parameters:
k = 16; % Number of sensors
m = 4; % Size of observation vectors b
n = 20; % Size of unknown vector x
reliable_sensors_list = [6 8 10 12 14]; % Number of consistent sensors
delta = 1e-6; % Concave approximation related constant
threshold = 1e-4; % Threshold to recover reliable sensors
methods = 4; % Methods being studied ( LS, l1, P1, P2(1) )

% save results
results = zeros(methods, length(reliable_sensors_list));
%preallocations
bi = zeros(m, 1, k);


for s_index = 1:length(reliable_sensors_list)
    
    s = reliable_sensors_list(s_index);
    
    % unknown vector is modeled as x0 ~ N(0, n^(-1/2)In)
    x0 = mvnrnd(zeros(1, n), n^(-0.5)*eye(n))';
    
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
    
    % Rearrange arrays and matrices
    b = bi(:);
    C = permute(Ai, [1 3 2]);
    A = reshape(C, [], size(Ai, 2), 1);
    
    reliable_sensors = [ones(1, s) zeros(1, k-s)];
    
    % LS method
    x_ls = ls_method(A, b, n);
    results_ls = sensor_validation(Ai, bi, x_ls, threshold, k, s);
    results(1, s_index) = results(1, s_index) + isequal(reliable_sensors, results_ls);
    
    % l1 method
    x_l1 = l1_method(A, b, n);
    results_l1 = sensor_validation(Ai, bi, x_l1, threshold, k, s);
    results(2, s_index) = results(2, s_index) + isequal(reliable_sensors, results_l1);
    
    % P1 method
    x_p1 = p1_method(Ai, bi, n, k);
    results_p1 = sensor_validation(Ai, bi, x_p1, threshold, k, s);
    results(3, s_index) = results(3, s_index) + isequal(reliable_sensors, results_p1);
    
    % P2(1) method
    x_p2_1 = p2_1_method(Ai, bi, n, k, x_p1, delta);
    results_p2_1 = sensor_validation(Ai, bi, x_p2_1, threshold, k, s);
    results(4, s_index) = results(4, s_index) + isequal(reliable_sensors, results_p2_1);
end