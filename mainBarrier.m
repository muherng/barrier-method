% Initialisation:
%A = [-rand(1)*5, 0.5; rand(1)*7, -2*rand(1); rand(1), rand(1); -1, -1];
%b = [1; 2; 3 ; 4];
%c = [1 ; 1]
%T = 1

Q = eye(2); %quadratic term
c = [0;0]; %linear term
A = eye(2) %constraint matrix
b = [0.5;-0.5];
x0 = [1;1]
T = 1;

fun = @(x) x'*Q*x + c'*x;
Aeq = zeros(1,size(Q,2));
beq = 0;
lb = [];
ub = [];
nonlcon = [];
options = optimset('Display', 'off','Algorithm','interior-point');
result_ip = fmincon(fun,x0,-A,-b,Aeq,beq,lb,ub,nonlcon,options);



% Newton's method parameters
%maxIter = 150 ;
maxIter = 150;
tol = 1e-10 ;
%tol = 0.15;
%alpha = 0.1 ;
alpha = 0.1;
beta = 0.5 ;


%% Interior points parameters
thresholdIP = 1e-5 %threshold for interior point method
%mu = 1.5
mu = 1.5;
%x0 = [1, 1]' %initial guess
T = 1;
dualityGap = zeros(1, maxIter) ; %track duality gap
T_values = zeros(1, maxIter) ; %track values of T parameter
N_Iterations = [] ; %used for the plotting
prev_i = 0 ; %used for the plotting
ord = [] ; %used for the plotting


%% Solve the program using Interior Point method with Newton + backtracking at each step
for j = 1:maxIter
    T_values(1,j) = T ;
    [x,i,lambda, f, dualf, dualArg] =Newton2(Q, A, b, c, T, x0, maxIter, tol, alpha, beta) ;
    dualityGap(:,j) = f - dualf ;
    N_Iterations = [N_Iterations, (1:i) + prev_i] ; %for plotting
    prev_i = prev_i + i ; %for plotting
    ord = [ord, repmat(dualityGap(:,j), 1, i)] ; %for plotting
    if (dualityGap(:,j) < thresholdIP)
        disp('Interior Point Threshold Achieved -- Exiting');
        break;
    else
        T = mu*T;
        x0 = x(:,i);
    end;
end;


%% Plot duality gap VS iterations number
% Figure1=figure(1);clf;
% set(Figure1,'defaulttextinterpreter','latex');
% stairs(N_Iterations, ord)
% xlabel('Iterations')
% ylabel('Duality Gap')
% title('Duality Gap depending on the number of iterations')
