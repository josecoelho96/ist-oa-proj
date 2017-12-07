% Least squares method, cvx aided
function [ x ] = ls_method( A, b, n )
    cvx_begin quiet
        variable x(n, 1)
        minimize( (b-A*x)'*(b-A*x) )
    cvx_end
end