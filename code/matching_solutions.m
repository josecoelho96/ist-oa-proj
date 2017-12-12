function [ L ] = matching_solutions( A, b, n, k, delta, x0, L0)

    cvx_begin quiet
        variable x(n, k) 
        variable L(n, 1)
        % define cost function
        f = 0;
        for i=1:k
            f = f + ((x(:,i)-L)'*(x(:,i)-L))/((x0(:,i)-L0)'*(x0(:,i)-L0));
        end
        minimize(f)
        subject to
        for i=1:k
            (b(:,:,i)-A(:,:,i)*x(:,i))'*(b(:,:,i)-A(:,:,i)*x(:,i)) <= delta;
        end    
    cvx_end
end