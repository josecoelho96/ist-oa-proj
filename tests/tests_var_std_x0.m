close all;
clear;

k = 16;
m = 4;
n = 20;

iter = 1000;
total = 100;

for j=1:total
    
    for i=1:iter
        % unknown vector is modeled as x0 ~ N(0, n^(-1/2)In)
        x0(:,i) = mvnrnd(zeros(1, n), n^(-0.5)*eye(n))';
        x0_prof(:,i) = randn(n,1) * 1/sqrt(n);
    end
    var_x0(:, j) = var(x0, 0, 2);
    var_x0_prof(:, j) = var(x0_prof, 0, 2);
end

% figure;
% hold on;
% for i=1:n
%     plot(var_x0(i,:))
% end

mean_var_x0 = mean(var_x0, 2);
mean_var_x0_prof = mean(var_x0_prof, 2);

fill = ones(total, 1);
figure;
hold on;
for i=1:n
    plot(fill*mean_var_x0(i))
end


figure;
hold on;
for i=1:n
    plot(fill*mean_var_x0_prof(i))
end