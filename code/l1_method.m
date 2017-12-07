% Least squares method (norm 1), cvx aided
function [ x ] = l1_method( A, b, n )
    cvx_begin quiet
        variable x(n, 1)
        minimize( norm(b-A*x, 1) )
    cvx_end 
end