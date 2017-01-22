function [V] = Newton_DD(X,Y,U)
%% Written by: Kapil Ahuja
%            (2014MT60663)
% The function implements the divided difference method to find an
% interpolating polynomial. Coefficients are generated in reverse order as
% the coeff at index i+1 is needed in calculating updated value at index i.
% 
% The multiplying matrix has terms of the form
% 1,(u(1)-x1),(u(1)-x2)*(u(1)-x1)....
% 1,(u(2)-x1),(u(2)-x2)*(u(2)-x1)....
% . . . . . . . . 
% . . . . . . . .
% So the values V are found by V = Matrix * coeff.
%%
n = length(X);
coeff = Y;
values = zeros(n);
for i = 1:n
    values(i,1) = Y(i);
end
for i = 2:n
    for j = 2:i
        values(i,j) = (values(i,j-1) - values(i-1,j-1))/(X(i) - X(i-j+1));
    end
end
coeff = diag(values)
%for i = 1:n-1
%    for j = n:-1:i+1
        %difference at i + 1, is used at i, hence the reverse order.
%        coeff(j) = (coeff(j) - coeff(j-1))/(X(j) - X(j-i));
%    end
%end
m = length(U);
VD = zeros(m,n);
%% Generate Van derMonde and find values at U.
for i = 1:m
    VD(i,1) = 1;
    for j = 2:n                                                 %generate Van derMonde matrix
        VD(i,j) = VD(i,j-1)*(U(i) - X(j-1));
    end
end
V = VD*coeff;