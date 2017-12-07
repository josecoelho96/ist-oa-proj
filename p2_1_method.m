% concave approximation method, cvx aided
function [ x ] = p2_1_method( A, b, n, k, x0, delta )
    for i=1:k
       w(i) = 1/(norm(b(:,:,i)-A(:,:,i)*x0) + delta);
    end
    cvx_begin quiet
        variable x(n, 1)
        % define cost function
        f = 0;
        for i=1:k
            f = f + w(i)*norm(b(:,:,i)-A(:,:,i)*x);
        end
        minimize(f)
    cvx_end
end