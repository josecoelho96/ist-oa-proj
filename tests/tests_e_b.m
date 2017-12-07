% generate bi and check expected values
close all;
clear;

m=4;
n=20;
k=16;
total = 1000;

norm_bi_r_opt1 = zeros(total, 1);
norm_bi_r_opt2 = zeros(total, 1);
norm_bi_nr = zeros(total, 1);

var_x0_opt1 = zeros(total, 1);
var_x0_opt2 = zeros(total, 1);

for i=1:total
    % A ~ N(0,1)
    Ai = randn(m, n);
    
    % x0 ~ N(0, n^(-1/2))
    x0_opt1 = mvnrnd(zeros(1, n), n^(-1)*eye(n))';
    % desvio padrao = 1/sqrt(n) !!!
    % variancia = 1/n !!!
    x0_opt2 = randn(n,1) * 1/sqrt(n);
        
    var_x0_opt1(i) = var(x0_opt1);
    var_x0_opt2(i) = var(x0_opt2);
    
    % reliable sensor
    bi_r_opt1 = Ai*x0_opt1;
    bi_r_opt2 = Ai*x0_opt2;
    % unreliable sensor
    bi_nr = mvnrnd(zeros(1, m), eye(m))';
        
    norm_bi_r_opt1(i) = norm(bi_r_opt1)^2;
    norm_bi_r_opt2(i) = norm(bi_r_opt2)^2;
    norm_bi_nr(i) = norm(bi_nr)^2;    
end

% check if expected values are similar
fprintf('Expected value for:\nReliable (Option 1): %f\nReliable (Option 2): %f\nNon reliable: %f\n\n', ...
mean(norm_bi_r_opt1), mean(norm_bi_r_opt2), mean(norm_bi_nr));

% check if variance is similar
fprintf('Mean variance value for:\nReliable (Option 1): %f\nReliable (Option 2): %f\n', ...
mean(var_x0_opt1), mean(var_x0_opt2));