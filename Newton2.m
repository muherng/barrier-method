
function [x,i,lambda, f, dualf, dualArg] = Newton2(Q, A, b, c, T, x0, maxIter, threshold, alpha, beta, quadratic_tail, penetrate)
%compute initial parameters
[m,n] = size(A);
x = zeros(n, maxIter);
x(:,1) = x0;
dx = [];
hess = [];
g = [];

q = -10; %Parameter set in paper  
s0 = -1/q;
C = [[2*s0, 1, 0];[2, 0, 0];[0, 0, 1]];
r = [-1/s0;1/(s0^2);-log(s0)];
a = C\r;

%begin iteration loop
for i = 1:maxIter
    %gradient and hessian
    %f = T*c'*x(:,i) - sum(log(b-A*x(:,i)));
    if A*x(:,i) >= b
        quadratic_tail = 0
    end
    
    if quadratic_tail
        A_pen = A(1:penetrate,:);
        s_pen = A_pen*x(:,i);
        quad_log = s_pen >= s0*ones(penetrate,1);
        quad_log = [quad_log; ones(m-penetrate,1)];
        A_feasible = zeros(m,n);
        A_unfeasible = zeros(m,n);
        A_feasible(find(quad_log),:) = A(find(quad_log),:)
        A_unfeasible(find(~quad_log),:) = A(find(~quad_log),:);
        b_feasible = zeros(m,1);
        b_unfeasible = zeros(m,1);
        b_feasible(find(quad_log)) = b(find(quad_log));
        b_feasible(find(~quad_log)) = b(find(~quad_log));
        feasible = sum(log(A_feasible * x(:,i)));
        unfeasible = sum(log(A_unfeasible*x(:,i)));
        f = T*x(:,i)'*Q*x(:,i) + T*c'*x(:,i) - feasible - unfeasible;
        d_feasible = 1./(A_feasible*(x:,i) - b_feasible)
        grad = T*diag(diag(Q))*x(:,i) + T*Q*x(:,i) + T*c - A_feasible'*d_feasible + ones(m,1)'*2*a(1)*diag(A_unfeasible*x(:,i)*A_unfeasible + a(2)*A_unfeasible);
        Hess = T*diag(diag(Q)) + T*Q + A_feasible' *diag(d_feasible.^2)*A_feasible + 2*a(1)*(A_unfeasible'*A_unfeasible)  
    else 
        f = T*x(:,i)'*Q*x(:,i) + T*c'*x(:,i) - sum(log(A*x(:,i) - b));
        d = 1./(A*x(:,i) - b);
        grad = T*diag(diag(Q))*x(:,i) + T*Q*x(:,i) + T*c - A'*d;
        Hess = T*diag(diag(Q)) + T*Q + A' * diag(d.^2) * A;
    end
    %disp('debug variables')
    %disp(T);
    %disp(c);
    %disp(x);
    %disp(T*c'*x(:,i));
    %disp(sum(log(b-A*x(:,i))));
    %disp(f);
 
    
    %d = 1./(b - A * x(:,i));
    %grad = T*c + A' * d;
    %Hess = A' * diag(d.^2) * A;
    
    %d = 1./(A*x(:,i) - b);
    %grad = T*diag(diag(Q))*x(:,i) + T*Q*x(:,i) + T*c - A'*d;
    %Hess = T*diag(diag(Q)) + T*Q + A' * diag(d.^2) * A;
    hess = [hess Hess];
    g = [g grad];
    
    
    %Newton step and decrement
    Dx = - Hess \ grad;
    lambda = - grad' * Dx;
    dx = [dx Dx];
    
    %stopping criterion
    if (lambda / 2) < threshold || i == maxIter || ~isempty(find(isnan(Dx),1)) || ~isempty(find(isinf(Dx),1))
        if (lambda / 2) < threshold
            disp('threshold exit')
        elseif i == maxIter
            disp('iteration exit')
        elseif  ~isempty(find(isnan(Hess),1)) || ~isempty(find(isinf(Hess),1))
            disp('boundary collision exit')
        end
        dualf = f -m/T;
        dualArg = 1/T * d;
        break
    else
        %backtracking line search
        bound = b - A*x(:,i);
        t = 1;
        %arbitrary_stepsize = 1;
        adx = A*Dx;
        for j = 1:size(bound,1)
            if ~isnan(bound(j)/adx(j)) && ~isinf(bound(j)/adx(j))
                if bound(j)/adx(j) > 0
%                     if arbitrary_stepsize == 1
%                         t = min(t,bound(j)/adx(j));
%                         arbitrary_stepsize = 0;
%                     else
%                         t = min(t,bound(j)/adx(j));
%                     end
                    t = min(t,bound(j)/adx(j));
                end
            end
        end
        t = t*.999;
        while T*(x(:,i)+t*Dx)'*Q*(x(:,i)+t*Dx) + T*c'*(x(:,i)+t*Dx) - sum(log(A*(x(:,i)+ t*Dx) - b)) > (f+alpha*t*grad'*Dx)
            t = beta * t;
        end 
        %update
        x(:,i+1) = x(:,i) + t * Dx ;
        if isnan(x(1,i+1)) || isinf(x(1,i+1))
            disp(t)
        end
    end 
end 
end
