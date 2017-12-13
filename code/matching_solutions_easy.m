function [ L ] = matching_solutions_easy( A, b, n, k, delta)
% CVX problem solving
cvx_begin quiet
    variable x(n, k)
    variable L(n, 1)
    % define cost function
    for i=1:k
        % f(i) = (x(:,i)-L)'*(x(:,i)-L);
        f(i) = norm(x(:,i)-L);
    end
    minimize(sum(f))
    subject to
    for i=1:k
        norm(bi(:,:,i) - Ai(:,:,i)*x(:,i)) <= restriction_delta;
        % (bi(:,:,i)-Ai(:,:,i)*x(:,i))'*(bi(:,:,i)-Ai(:,:,i)*x(:,i)) <= restriction_delta;
    end
cvx_end
end