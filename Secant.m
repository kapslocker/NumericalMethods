function s1 = Secant(x0,x1,in1)
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
s1 = x0;
x = [x0];
y = [f1(x0)];
flag = 0;                       % to check if solution has been found.
epsilon = 1e-8;                 % to check if |grad(f)| is very low. 
maxit = 20;
tol = 1e-8;                     % The relative tolerance beyond which answer is printed.
flag2 = 0;                      % This is for checking grad(f)'s termination.
for i=1:maxit
    temp = f1(x1);
    temp_grad = (f1(x1) - f1(x0))/(x1 - x0);
    x = [x x1];
    y = [y temp];
    if(abs(temp_grad) < epsilon)
        disp 'The value of grad(f) is too small to divide, and will cause errors.';
        flag2 = 1;
        break;
    end
    if(abs(x1-x0) <= tol*abs(x0) )
        flag = 1;
        s1 = x1;
        break;
    end
    x0 = x1;
    x1 = x1 - temp/temp_grad;
end
%%Plotting solution now.
if(flag && ~isequal(length(x),0))
    X = sprintf('Converged after %d iterations',i);
    disp(X);
    l = min(x);
    u = max(x);
    t = [l-1:0.1:u+1];
    z = f1(t);
    plot(x,y,t,z,'.-');
    xlabel('x')
    ylabel('f(x)')   
end
if(~flag)
    if(~flag2)
        disp 'The function is diverging'
    end
    l = min(x);
    u = max(x);
    t = [l-1:0.5:u+1];
    plot(x,y,t,f1(t),'.-');
    xlabel('x')
    ylabel('f(x)') 
end