% generate bi and check expected values
close all;
clear;

m=4;
n=20;
k=16;

% SNR = 10*log10(sigma^(-2))
SNR = 5;
% sigma = 10^(-SNR/20)
sigma = 10^(-SNR/20);

total = 100000;

norm_bi_r = zeros(total, 1);
norm_bi_r_opt2 = zeros(total, 1);
norm_bi_nr = zeros(total, 1);
norm_bi_nr_opt2 = zeros(total, 1);

for i=1:total
    
    
    % A ~ N(0,1)
    Ai = randn(m, n);
    
    % x0 ~ N(0, n^(-1/2))
    x0 = mvnrnd(zeros(1, n), n^(-1)*eye(n))';
    
    % vi ~ N (0, sigma^2)
    vi = mvnrnd(zeros(1, m), (sigma^2)*eye(m))';
    vi_opt2 = randn(m, 1)*sigma;

    % reliable sensor
    bi_r = Ai*x0 + vi;
    bi_r_opt2 = Ai*x0 + vi_opt2;
    
    % unreliable sensor
    bi_nr = mvnrnd(zeros(1, m), (1+sigma^2)*eye(m))';
    bi_nr_opt2 = randn(m, 1) * sqrt(1+sigma^2);
    
    norm_bi_r(i) = norm(bi_r)^2;
    norm_bi_r_opt2(i) = norm(bi_r_opt2)^2;
    norm_bi_nr(i) = norm(bi_nr)^2;    
    norm_bi_nr_opt2(i) = norm(bi_nr_opt2)^2;    
end


% check if expected values are similar
fprintf('[Option 1] Expected value for:\nReliable: %f\nNon reliable: %f\n\n', ...
mean(norm_bi_r), mean(norm_bi_nr));

fprintf('[Option 2] Expected value for:\nReliable: %f\nNon reliable: %f\n\n', ...
mean(norm_bi_r_opt2), mean(norm_bi_nr_opt2));