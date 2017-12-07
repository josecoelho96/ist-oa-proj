close all;
clear;

n=20;

% representation
% N(mu, sigma^2)
% mu : media
% sigma^2 : variancia
% sigma : desvio padrao 
% desvio padrao = sqrt(variancia)

% N(media, variancia)
% variancia = n^(-1/2)

number_tests = 100;

var_desired = 1/sqrt(n);
std_desired = sqrt(var_desired);

std_des = ones(1, number_tests) * std_desired;
var_des = ones(1, number_tests) * var_desired;

std_sug = zeros(1, number_tests);
std_op1 = zeros(1, number_tests);
std_op2 = zeros(1, number_tests);

var_sug = zeros(1, number_tests);
var_op1 = zeros(1, number_tests);
var_op2 = zeros(1, number_tests);

for i=1:number_tests    
    %sugestao prof
    x0_i = randn(n,1) * 1/sqrt(n);
    %opcao 1
    x0_t1 = mvnrnd(zeros(1, n), n^(-0.5)*eye(n));
    %opcao 2
    x0_t2 = mvnrnd(zeros(1, n), n^(-1)*eye(n));
    
    std_sug(i) = std(x0_i);
    std_op1(i) = std(x0_t1);
    std_op2(i) = std(x0_t2);

    var_sug(i) = var(x0_i);
    var_op1(i) = var(x0_t1);
    var_op2(i) = var(x0_t2);
end

mean_std_sug = mean(std_sug);
mean_std_op1 = mean(std_op1);
mean_std_op2 = mean(std_op2);

mean_var_sug = mean(var_sug);
mean_var_op1 = mean(var_op1);
mean_var_op2 = mean(var_op2);

m_std_sug = ones(1, number_tests)*mean_std_sug;
m_std_op1 = ones(1, number_tests)*mean_std_op1;
m_std_op2 = ones(1, number_tests)*mean_std_op2;

m_var_sug = ones(1, number_tests)*mean_var_sug;
m_var_op1 = ones(1, number_tests)*mean_var_op1;
m_var_op2 = ones(1, number_tests)*mean_var_op2;

figure;
hold on;
plot(std_sug, 'r');
plot(std_op1, 'b');
plot(std_op2, 'g');
title('Standard deviation obtained');
legend('Prof sug.', 'Op1', 'Op2');
ylim([0 0.8]);

figure;
hold on;
plot(std_des, 'k');
plot(m_std_sug, 'r');
plot(m_std_op1, 'b');
plot(m_std_op2, 'g');
title('Average standard deviation');
legend('Paper', 'Prof sug.', 'Op1', 'Op2');
ylim([0 0.8]);



figure;
hold on;
plot(std_sug, 'r');
plot(std_op1, 'b');
plot(std_op2, 'g');
title('Standard deviation obtained');
legend('Prof sug.', 'Op1', 'Op2');
ylim([0 0.8]);


figure;
hold on;
plot(var_des, 'k.');
plot(m_var_sug, 'r');
plot(m_var_op1, 'b');
plot(m_var_op2, 'g');
title('Average variance');
legend('Paper', 'Prof sug.', 'Op1', 'Op2');
ylim([0 0.5]);