% convex relaxation method, cvx aided
function [ x ] = p1_method( A, b, n, k )
    cvx_begin quiet
        variable x(n, 1)
        % define cost function
        f = 0;
        for i=1:k
            f = f + norm(b(:, :, i) - A(:, :, i)*x);
        end
        minimize(f)
    cvx_end
end