% Matching solutions method
function [ x, L ] = matching_solutions_miter( A, b, n, k, delta, x0, L0)

    cvx_begin quiet
        variable x(n, k) 
        variable L(n, 1)
        % define cost function
        for i=1:k
            f(i) =  norm(x(:,i)-L) / norm(x0(:,i)-L0);
        end

        minimize(sum(f))

        subject to
        for i=1:k
            (b(:,:,i)-A(:,:,i)*x(:,i))'*(b(:,:,i)-A(:,:,i)*x(:,i)) <= delta;
        end
    cvx_end
end