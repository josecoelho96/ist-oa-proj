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
s = 14; % Number of consistent sensors
delta = 1e-6; % Concave approximation related constant
threshold = 1e-4; % Threshold to recover reliable sensors

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

% LS method
x_ls = ls_method(A, b, n);
results_ls = sensor_validation(Ai, bi, x_ls, threshold, k, s);

% l1 method
x_l1 = l1_method(A, b, n);
results_l1 = sensor_validation(Ai, bi, x_l1, threshold, k, s);

% P1 method
x_p1 = p1_method(Ai, bi, n, k);
results_p1 = sensor_validation(Ai, bi, x_p1, threshold, k, s);

% P2(1) method 
x_p2_1 = p2_1_method(Ai, bi, n, k, x_p1, delta);
results_p2_1 = sensor_validation(Ai, bi, x_p2_1, threshold, k, s);