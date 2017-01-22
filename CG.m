function [x] = CG(A,b,M,x0)
%Written by: Kapil Ahuja
%            (2014MT60663)
% Since a preconditioner is given, the norm used is the M norm, which is
% given by z'*r (instead of r'*r). Here z is the solution of Mz = r. Since
% CG method always gives a solution in under n iterations, n being the
% dimensions of the matrix, we loop for n interations. The direction is
% p, at each step and is updated as p = z + (z'*r/prev(z'*r)) * p ( The 
% rest is as  derived in lecture slides). 
%%
s = size(A);
n = s(1);
r = b - A*x0;
z = M\r;
er = (z')*r;
e0 = er;
p = z;
x = x0;
tol = 1e-6;
for k = 1:n
    s = A*p;
    alpha = er/((p')*s);
    x = x + alpha*p;                                                        % update x in the new iteration.
    r = r - alpha*s;
    z = M\r;
    temp = er;
    er = (z')*r;                                                            % update error.
    if(er<tol*e0)                                                           % check for solution.
        break
    end
    beta = er/temp;
    p = z + beta*p;                                                         %update direction.
end