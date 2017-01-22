function [V] = cubic_spline(X,Y,U)
%Written by: Kapil Ahuja
%            (2014MT60663)
% The interpolating function in interval [x i,x i+1] is given by
% s(x)=(1-t)*f(xi) + t*f(x i+1) + t*(1-t)*(a*(1-t) + b*t)
% where : t = x - xi /(x i+1 - xi)
% a = ki*(xi+1 - xi) - (f(x i+1) - f(x i))
% b = -k i+1 * (x i+1 - xi ) + (f(x i+1) - f(x i))
% ki = f_der(xi), k i+1 = f_der(x i+1)
% We solve for ki's first, which then determine ai's and bi's.
% ki's are obtained by solving a linear system. The matrix A is tridiagonal
% Later ai's, bi's are calculated and are evaluated by applying binary search. 
%%
n = length(X);
[X,I] = sort(X);
Y = Y(I);
A = zeros(n,n-2);
b = zeros(1,n-2);
%% storing 4n-2 conditions first.
for i = 1:n-2
    A(i,i) = X(i+2) - X(i+1);
    A(i+1,i) = 2 *( (X(i+2) - X(i)) + (X(i)) );
    A(i+2,i) = X(i+1) - X(i);
    d1 = X(i+2) -X(i+1);
    d2 = X(i+1) - X(i);
    d3 = Y(i+2) - Y(i+1);
    d4 = Y(i+1) - Y(i);
    b(i) = 3 * (d1*d4/d2 + d3*d2/d1 );
end
%% Adding the natural spline conditions and resizing the matrix.
A = [zeros(n,1) A zeros(n,1)];
A(1,1) = 2;
A(2,1) = 1;                                                 
A (n-1,n) = 1;
A(n,n) = 2;

b = [0,b,0];    
b(1) = 3* (Y(2) - Y(1) )/ ( X(2) - X(1) );
b(n) = 3* (Y(n) - Y(n-1) )/ ( X(n) - X(n-1) );

%% Now solve for k.
K = b/inv(A);

%% Find ai's and bi's correspoding to given ai's.
a = zeros(n,1);
b = zeros(n,1);
for i =1:n-1
    a(i) = K(i)*(X(i+1) - X(i)) - (Y(i+1) - Y(i));
    b(i) = -K(i+1)*(X(i+1) - X(i))  +(Y(i+1) - Y(i));
end
%% Apply binary search to find the correct interval and evaluate interpolated values. 
m = length(U);
n = length(X);
V = zeros(m,1);
for i = 1:m
    low = 1;
    high = n;
    pos = -1;
    while(low<high)
        mid = floor((low + high)/2);
        if(X(mid)<= U(i) && X(mid + 1)>= U(i) )             % This is the right position for 
            pos = mid;                                      % U(i) and we need the interpolating 
            break;                                          % polynomial in this interval.
        end
        if(X(mid+1) < U(i))
            low = mid + 1;
        elseif(X(mid)>U(i))
            high = mid;
        end
    end
    t = (U(i) - X(pos)) / (X(pos+1) - X(pos));
    V(i) = (1 -t)*Y(pos) + t*Y(pos+1) + t*(1-t)*(a(pos)*(1-t) + b(pos)*t) ;
end
%% plot the found points.
plot(U,V,X,Y,'.-');
xlabel('x')
ylabel('f(x)')
