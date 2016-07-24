% Custom Interior Point Solver
% Quadratic barrier for unfeasible initial guess
% Modified Newton Step 
% Todorov 2011 -- A Convex Smooth and Invertible Contact Model for Trajectory Optimization 

%Problem Statement
%Minimize x'Qx + c'x 
%Subject to Ax >= b
%           lb <=  x <= ub

function [x_min,f,feasible] = InteriorPoints(Q, A, b, c, ub, lb, x0)

%% Newton's method parameters
num_params = size(A,2);
%[penetrate,~] = size(b);
A_full = full([A; eye(num_params); -eye(num_params)]);
b_full = [b; lb; -ub];

maxIter = 3;  %was 3
newtonIter = 5; %was 4
alpha = 0.01;
beta = 0.5 ;

%% Interior points parameters
T = 1; %inverse of log barrier constant
thresholdIP = 1e-5; %threshold for interior point method
mu = 100; %used to be 100  %schedule of T (consider revising this)
re = 1.0;
dualityGap = zeros(1, maxIter) ; %track duality gap
T_values = zeros(1, maxIter) ; %track values of T parameter
quadratic_tail = 0;

%% Find Feasible Initial Guess
% if ~isempty( find(A*x0 <= b) )
%    quadratic_tail = 1; 
% end

%% Solve the program using Interior Point method with Newton + backtracking at each step
for j = 1:maxIter
    T_values(1,j) = T ;
    [x,i,lambda, f, dualf, dualArg, feasible] =Newton2(Q, A_full, b_full, c, T, x0, maxIter, thresholdIP, alpha, beta, newtonIter);
    dualityGap(:,j) = f - dualf ;
    
    %I'm going to comment termination condition out, but remember to think
    %more about an appropriate termiantion condition
%     if (dualityGap(:,j) < thresholdIP)
%         disp('Are you kidding me?');
%         break;
%     else
%         T = mu*T;
%         x0 = x(:,i);
%     end;
    if feasible
        T = mu*T;
    else
        T = (1/re)*T;
    end
    x0 = x(:,i);
    
end;
x_min = x(:,i);

end
