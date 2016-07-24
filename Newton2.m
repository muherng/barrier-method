
function [x,i,lambda, f, dualf, dualArg, feasible] = Newton2(Q, A, b, c, T, x0, maxIter, threshold, alpha, beta, newtonIter)
%compute initial parameters
[m,n] = size(A);
x = zeros(n, maxIter);
x(:,1) = x0;
dx = [];
hess = [];
g = [];

%quadratic coefficients
q = -10; %Parameter set in paper  
s0 = -1/q;
C = [[2,0,0]; [2*s0, 1, 0];[0, 0, 1]]; %quadratic, linear, and constant term
r = [1/(s0^2);-1/s0;-log(s0)];
a = C\r;

quad = 1; %need quadratic tail
%begin iteration loop
for i = 1:newtonIter
    %gradient and hessian
    %f = T*c'*x(:,i) - sum(log(b-A*x(:,i)));
    if isempty(find(A*x(:,i) <= b))
        quad = 0;
    else
        quad = 1;  %but I don't think this needs to be done.  
    end
    
    if quad
        s = A*x(:,i) - b;
        %disp(s);
        quad_log = ~(s <= s0);
        %binary vector, 1 for log evaluation, 0 for quadratic
        Af = A(find(quad_log),:);
        Au = A(find(~quad_log),:);
        bf = b(find(quad_log));
        bu = b(find(~quad_log));
        feasible = sum(log(Af * x(:,i) - bf));
        unfeasible = sum(horzcat((Au*x(:,i)-bu).^2,Au*x(:,i) - bu,ones(size(find(~quad_log))))*a);
        %f = T*x(:,i)'*Q*x(:,i) + T*c'*x(:,i) - feasible + unfeasible;
        f = - feasible + unfeasible;
        d = 1./(Af*x(:,i) - bf); 
        %grad = T*(Q + Q')*x(:,i) + T*c - Af'*d + (2*a(1)*diag(Au*x(:,i)-bu)*Au + a(2)*Au)'*ones(size(find(~quad_log)));
        %Hess = T*(Q + Q') + Af'*diag(d.^2)*Af + 2*a(1)*(Au'*Au);
        grad = - Af'*d + (2*a(1)*diag(Au*x(:,i)-bu)*Au + a(2)*Au)'*ones(size(find(~quad_log)));
        Hess = Af'*diag(d.^2)*Af + 2*a(1)*(Au'*Au);
        %f = sum(horzcat((A*x(:,i)-b).^2,A*x(:,i) - b,ones(m,1))*a);
        %grad = (2*a(1)*diag(A*x(:,i)-b)*A + a(2)*A)'*ones(m,1);
        %Hess = 2*a(1)*(A'*A);
    else 
        f = T*x(:,i)'*Q*x(:,i) + T*c'*x(:,i) - sum(log(A*x(:,i) - b));
        d = 1./(A*x(:,i) - b);
        %grad = T*diag(diag(Q))*x(:,i) + T*Q*x(:,i) + T*c - A'*d;
        grad = T*(Q + Q')*x(:,i) + T*c - A'*d;
        Hess = T*(Q + Q') + A' * diag(d.^2) * A;
        %Hess = T*diag(diag(Q)) + T*Q + A' * diag(d.^2) * A;
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
    
%     disp('here is x');
%     disp(x(:,i));
%     disp('here is grad');
%     disp(grad);
%     disp('here is Hess');
%     disp(Hess);
%     disp('descent direction');
%     disp(Dx);
    
    %stopping criterion
    %if (lambda / 2) < threshold || i == newtonIter || ~isempty(find(isnan(Dx),1)) || ~isempty(find(isinf(Dx),1))
    if (lambda / abs(f)) < threshold || i == newtonIter || ~isempty(find(isnan(Dx),1)) || ~isempty(find(isinf(Dx),1))
        if (lambda / abs(f)) < threshold
            if quad
                disp('threshold exit');
                disp(i);
%                  disp('lambda');
%                  disp(lambda);
%                  disp('f');
%                  disp(f);
%                 disp('Dx');
%                 disp(Dx);
%                 disp('grad');
%                 disp(grad);
%                 disp('Hess');
%                 disp(Hess);
            end
            %disp(lambda);
        elseif i == maxIter
            if quad
                disp('iteration exit')
            end
        elseif  ~isempty(find(isnan(Hess),1)) || ~isempty(find(isinf(Hess),1))
            disp('boundary collision exit')
        end
        if quad == 1
            feasible = 0;
%             disp('unfeasible solution');
%             disp('right');
        else 
            feasible = 1;
%             disp('feasible solution');
        end
        %dualf = f -m/T;
        %dualArg = 1/T * d;
        dualf = 0;
        dualArg = 0;
        break
    else
        %backtracking line search
        t = 1;
        if ~quad
            bound = b - A*x(:,i);
            adx = A*Dx;
            for j = 1:size(bound,1)
                if ~isnan(bound(j)/adx(j)) && ~isinf(bound(j)/adx(j))
                    if bound(j)/adx(j) > 0
                        t = min(t,bound(j)/adx(j));
                    end
                end
            end
            t = t*.999;
        end
        backtrack = 1;
        %while T*(x(:,i)+t*Dx)'*Q*(x(:,i)+t*Dx) + T*c'*(x(:,i)+t*Dx) - sum(log(A*(x(:,i)+ t*Dx) - b)) > (f+alpha*t*grad'*Dx)
        fnew_list = [];
        f_list = [];
        while backtrack
            %here is what's happening, we can potentially descend quite a
            %bit, that's what the derivative projection suggests but
            %evaluation is not as good.  this is impossible, at some point,
            %unless the derivative is very small or if the derivative is
            %wrong
            
%              disp('stuck');
%               disp(quad);
%               disp('lambda');
%               disp(lambda);
%               disp('Dx');
%               disp(Dx);
%               disp('tgradDx');
%               disp(t*grad'*Dx);
%               disp('t');
%               disp(t)
              
              
%             disp(quad);
%             disp('here is x');
%             disp(x(:,i));
%             disp('descent direction');
%             disp(Dx);
%             disp('here is f');
%             disp(f);
%             disp('here is t')
%             disp(t);
            x_new = x(:,i) + t*Dx;
            if quad
                 feasible = sum(log(Af * x_new - bf));
                 unfeasible = sum(horzcat((Au*x_new-bu).^2,Au*x_new - bu,ones(size(find(~quad_log))))*a);
%                 %f_new = T*x_new'*Q*x_new + T*c'*x_new - feasible + unfeasible;
                 f_new = - feasible + unfeasible;
                  %f_new = sum(horzcat((A*x_new-b).^2,A*x_new - b,ones(m,1))*a);
            else
                f_new = T*x_new'*Q*x_new + T*c'*x_new - sum(log(A*x_new - b));
            end
            %if f_new < (f+alpha*t*grad'*Dx)
            if f_new < f + alpha*t*grad'*Dx
                backtrack = 0;
            else
                t = beta*t;
            end
%             f_list = [f_list f+alpha*t*grad'*Dx];
%             fnew_list = [fnew_list f_new];
%             disp(fnew_list(1:min(10,size(fnew_list,2))));
%             disp(f_list(1:min(10,size(f_list,2))));
            %disp(f_new);
            %disp(f )
           
        end
        %disp('unstuck');
        %update
        x(:,i+1) = x_new;
        if isnan(x(1,i+1)) || isinf(x(1,i+1))
            disp(t)
        end
    end 
end 
end
