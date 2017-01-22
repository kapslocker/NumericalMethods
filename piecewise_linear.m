function [V] = piecewise_linear(X,Y,U)
% Written by: Kapil Ahuja
%            (2014MT60663)
% Given any input x, first we find the interval x lies in. 
% To do this, binar search is used below. Since it has been guaranteed that
% the values lie in the range X is given, the condition is not handled
% later on. In an interval [xi , x i+1], the interpolating polynomial is the
% straight line connecting (xi,yi) and (x i+1 , y i+1).
%%
[X,I] = sort(X);
Y = Y(I);
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
    V(i) = Y(pos) + ( (Y(pos+1) - Y(pos)) * (U(i) - X(pos)))/(X(pos+1) - X(pos) ) ;
end