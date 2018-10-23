% this solves the second small problem
%
% minimize_X    (   f(X)  +  0.5*rho*(X - N)^2)
%
% where X = [A, B] and f(X) =  0 or infty depending on whehter y_i x_i' B >= 1 - A_i or not
%
%
%
function AB =  solve_small_problem_two(N_A, N_B,i,x,y)

%     c = pos(y(i)*(N_B'*x(:,i) + N_A(i) - 1)/(1 + x(:,i)'*x(:,i)));
%     
%     B = N_B - c*y(i)*x(:,i);
%     A = N_A - c;


    p = length(N_B);
    
    if (            y(i)*N_B'*x(:,i) >= 1 - N_A(i)     )
        A = N_A;
        B = N_B;
    else
        A = N_A;
        
        B = (eye(p)  -  (x(:,i)*x(:,i)')/(1 + x(:,i)'*x(:,i)) )*  (    N_B + (1 -  N_A(i))*y(i)*x(:,i)   );
        A(i) = 1 - y(i)*B'*x(:,i);
        
    end
    
    AB = [A;B];
    
end