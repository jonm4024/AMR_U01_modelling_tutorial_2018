%% Demo by Jose Bento @ BC 

%% create some random problem

n = 10;
p = 4000;
y = sign(randn(n,1)); 
x = randn(p,n);
x = [x;ones(1,n)];


%% solve using existing solver

cvx_begin

    variable gap(n,1);
    variable hyperplane(p+1);

    minimize    (       sum(  gap) +     0.5* (norm(hyperplane)^2)                  )

    subject to 
        
        for i =1:n
           y(i)*hyperplane'*x(:,i) >= 1 - gap(i);
        end
        gap >= 0;

cvx_end

%% solve using ADMM
%
%
%  our problem is    minimize   f_1 (X) + f_2 (X)  + ...  + f_{n+1} (X)  + f_{n+2} (X)
%
%  where X = [A, B], which has dimensions n  + p+1  and  
%  f_1(X) =  sum(A) + 0.5* (norm(B)^2)
%  f_{i+1}(X) =  0 or infty depending on whehter y_i x_i' B >= 1 - A_i or not
%  f_{n+2}(X) =  0 or infty depending on whehter A>=0 or not

rho = 1;
alpha = 1;

% initialize the variables
X1  = randn(n+p+1,1);
X2  = randn(n+p+1,n); %notice that X2 is actually n variables
X3  = randn(n+p+1,1);

Z = randn(n+p+1,1);

L1  = randn(n+p+1,1);
L2  = randn(n+p+1,n); %notice that L2 is actually n variables
L3  = randn(n+p+1,1);

% go over the ADMM iterations
% uncomment the next line and use parfor to exploit parallelism
% delete(gcp('nocreate')); parpool(2); 
for t = 1:2000
    disp(t)
    % do the X updates
    X1 = solve_small_problem_one(Z(1:n) - L1(1:n), Z(n+1:end) - L1(n+1:end), rho);

    % replace for by parfor for parallelism
    for i =1:n
        X2(:,i) = solve_small_problem_two(Z(1:n) - L2(1:n,i), Z(n+1:end) - L2(n+1:end,i), i,x,y);
    end
    X3 =  solve_small_problem_three(Z(1:n) - L3(1:n), Z(n+1:end) - L3(n+1:end));

    % do the Z updates
    Z = (L1+X1 + sum(L2+X2,2) + L3+X3) / (n + 2);

    % do the U updates
    L1 = L1 + alpha*(X1 - Z);
    % replace for by parfor for parallelism
    for i = 1:n
        L2(:,i) = L2(:,i) + alpha*(X2(:,i) - Z);
    end
    L3 = L3 + alpha*(X3 - Z);
end

max(abs(Z - [gap;hyperplane]))