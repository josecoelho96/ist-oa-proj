function [ L ] = matching_solutions_easy( A, b, n, k, delta)
    cvx_begin quiet
        variable x(n, k)
        variable L(n, 1)
        % define cost function
        f = 0;
        for i=1:k
            f = f + (x(:,i)-L)'*(x(:,i)-L);
        end
        minimize(f)
        subject to
        for i=1:k
            (b(:,:,i)-A(:,:,i)*x(:,i))'*(b(:,:,i)-A(:,:,i)*x(:,i)) <= delta;
        end    
    cvx_end
end