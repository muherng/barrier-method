% Custom Interior Point Solver
% Quadratic barrier for unfeasible initial guess
% Modified Newton Step 
% Todorov 2011 -- A Convex Smooth and Invertible Contact Model for Trajectory Optimization 

%Problem Statement
%Minimize x'Qx + c'x 
%Subject to Ax >= b
%           lb <=  x <= ub

function [x_min,f] = InteriorPoints(Q, A, b, c, ub, lb, x0)

%% Newton's method parameters
num_params = size(A,2);
[penetrate,~] = size(b);
A_full = full([A; eye(num_params); -eye(num_params)]);
b_full = [b; lb*ones(num_params,1); -ub*ones(num_params,1)];

maxIter = 150;
tol = 1e-10 ;
alpha = 0.1;
beta = 0.5 ;

%% Interior points parameters
T = 1; %inverse of log barrier constant
thresholdIP = 1e-5; %threshold for interior point method
mu = 1.5; %schedule of T (consider revising this)
dualityGap = zeros(1, maxIter) ; %track duality gap
T_values = zeros(1, maxIter) ; %track values of T parameter
quadratic_tail = 0;

%% Find Feasible Initial Guess
if A*x <= b
   quadratic_tail = 1; 
end

%% Solve the program using Interior Point method with Newton + backtracking at each step
for j = 1:maxIter
    T_values(1,j) = T ;
    [x,i,lambda, f, dualf, dualArg] =Newton2(Q, A_full, b_full, c, T, x0, maxIter, tol, alpha, beta, quadratic_tail,penetrate);
    dualityGap(:,j) = f - dualf ;
    if (dualityGap(:,j) < thresholdIP)
        %disp('Interior Point Threshold Achieved -- Exiting');
        disp('interior objective');
        disp(x(:,i)'*Q*x(:,i) + c'*x(:,i));
        break;
    else
        T = mu*T;
        x0 = x(:,i);
    end;
end;
x_min = x(:,i);

end
