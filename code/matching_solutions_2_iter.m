function [ L ] = matching_solutions_2_iter( A, b, n, k, delta, x0, L0)

%     cvx_begin quiet
%         variable x(n, k) 
%         variable L(n, 1)
%         % define cost function
%         f = 0;
%         for i=1:k
%             f = f + ((x(:,i)-L)'*(x(:,i)-L))/((x0(:,i)-L0)'*(x0(:,i)-L0));
%         end
%         minimize(f)
%         subject to
%         for i=1:k
%             (b(:,:,i)-A(:,:,i)*x(:,i))'*(b(:,:,i)-A(:,:,i)*x(:,i)) <= delta;
%         end    
%     cvx_end
    
    cvx_begin quiet
        variable x(n, k) 
        variable L(n, 1)
        % define cost function
        for i=1:k
            % f(i) = ((x(:,i)-L)'*(x(:,i)-L))/((x0(:,i)-L0)'*(x0(:,i)-L0));
            f(i) =  norm(x(:,i)-L) / norm((x0(:,i)-L0));
        end

        minimize(sum(f))

        subject to
        for i=1:k
            % norm(bi(:,:,i) - Ai(:,:,i)*x(:,i)) <= restriction_delta;
            % norm(bi(:,:,i) - Ai(:,:,i)*x(:,i)) <= sqrt(restriction_delta);
            (b(:,:,i)-A(:,:,i)*x(:,i))'*(b(:,:,i)-A(:,:,i)*x(:,i)) <= delta;
        end
    cvx_end
    
    x1 = x;
    L1 = L;
    
    cvx_begin quiet
        variable x(n, k) 
        variable L(n, 1)
        % define cost function
        for i=1:k
            % f(i) = ((x(:,i)-L)'*(x(:,i)-L))/((x0(:,i)-L0)'*(x0(:,i)-L0));
            f(i) =  norm(x(:,i)-L) / norm((x1(:,i)-L1));
        end

        minimize(sum(f))

        subject to
        for i=1:k
            % norm(bi(:,:,i) - Ai(:,:,i)*x(:,i)) <= restriction_delta;
            % norm(bi(:,:,i) - Ai(:,:,i)*x(:,i)) <= sqrt(restriction_delta);
            (b(:,:,i)-A(:,:,i)*x(:,i))'*(b(:,:,i)-A(:,:,i)*x(:,i)) <= delta;
        end
    cvx_end
    
    
end