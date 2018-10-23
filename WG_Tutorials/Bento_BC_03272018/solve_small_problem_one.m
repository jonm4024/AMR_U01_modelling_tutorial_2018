% this solves the first small problem
%
% minimize_X    (   f(X)  +  0.5*rho*(X - N)^2)
%
% where X = [A, B] and f(X) =  sum(A) + 0.5* (norm(B)^2)
%
%
%

function AB =  solve_small_problem_one(N_A, N_B, rho)

    A = N_A - (1/rho);
    B = 0.5*N_B;

    AB = [A; B];
    
end