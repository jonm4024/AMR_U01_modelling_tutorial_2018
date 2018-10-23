% this solves the second third problem
%
% minimize_X    (   f(X)  +  0.5*rho*(X - N)^2)
%
% where X = [A, B] and f(X) =  0 or infty depending on whehter A>=0 or not
%
%
%

function AB =  solve_small_problem_three(N_A, N_B)

    B = N_B;
    A = pos(N_A);
    
    AB = [A;B];
    
end