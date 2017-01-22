function sol = Newton(x0,in1,in2)
%% Written by: Kapil Ahuja
%            (2014MT60663)
% The function return the solution obtained by Newton's method on an
% initial guess x0 and the function defined in f.m and whose derivative is
% given in f_grad.m
% The stopping criteria is when relative error is less than tolerance.
% We do not divide the f'(x), when it gets too low(less than epsilon).
% The code runs for a maximum of 20 iterations, beyond which it is
% reasonable to assume that the function is diverging.
%
%%
f1 = inline(in1);
f1_der = inline(in2);
sol = x0;
x = [];
y = [];
flag = 0;
epsilon = 1e-8;
maxit = 20;
tol = 1e-5;
flag2 = 0;                      % This is for checking grad(f)'s termination.
for i=1:maxit
    temp = f1(x0);
    temp_grad = f1_der(x0);
    x = [x x0];
    y = [y temp];
    if(abs(temp_grad) < epsilon)
        disp 'The value of grad(f) is too small to divide, and will cause errors.';
        flag2 = 1;
        break;
    end
    x1 = x0 - temp/temp_grad;
    if(abs(x1-x0) <= tol*abs(x0))
        flag = 1;
        sol = x1;
        break;
    end
    x0 = x1;
end
%% print num_iterations and plot if solution found.
if(flag && ~isequal(length(x),0))
    X = sprintf('Converged after %d iterations',i);
    disp(X);
    l = x(end);
    u = x(1);
    t = [l-0.1:0.1:u+0.1];
    plot(x,y,t,f1(t));
    xlabel('x')
    ylabel('f(x)')
end
%% This is when solution is not found or the iterations have diverged.
if(~flag)
    if(~flag2)
        disp 'The function is diverging'
    end
    l = min(x);
    u = max(x);
    t = [l-1:0.1:u+1];
    x
    plot(x,y,t,f1(t),'.-');
    xlabel('x')
    ylabel('f(x)') 
end