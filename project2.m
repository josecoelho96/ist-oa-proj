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
s = 6; % Number of consistent sensors

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
