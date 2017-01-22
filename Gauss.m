function x = Gauss(A,b)
%
% Written by: Kapil Ahuja
%            (2014MT60663)
% Performs Gauss elimination on an Augmented Matrix (A|b). Here A has
% dimensions m*n and b has dimensions m*1. The pivoting used here is that
% the element with the maximum absolute value in the k'th iteration is
% brought to the diagonal position. Invertibility is checked mid way, which
% fails when at any iteration, there is no element with a non zero absolute
% value, to replace the diagonal element(which is also zero).
%
%
%%
% Augment the matrices, pivot and convert to Upper matrix.
A = [A b];
s = size(A);
m = s(1);
n = m+1;
for k = 1:m
    i_max = k;
    for i = k:m                                   %pivot
        if( abs(A(i,k)) > abs(A(i_max,k)) )
            i_max = i;
        end
    end
    if(A(i_max,k) == 0)                          % check for invertibility
        error('Matrix is non-invertible.');
    end
    A([k,i_max],:) = A([i_max,k],:);             % swap rows k and i_max of augmented matrix.
    for i = k+1:m
        ratio = A(i,k)/A(k,k);
        for j = k+1:n
            A(i,j) = A(i,j) - ratio*A(k,j);
        end
        A(i,k)=0;
    end
end
%% Solve using backward substitution now.

b = A(:,end);                                    % extracting b(modified) from the Augmented Matrix
A = A(:,1:end-1);
s = size(A);n = s(2);
for k = n:-1:1
    x(k) = b(k);
    for j = k+1:n
        x(k) = x(k) - A(k,j)*x(j);
    end
    x(k) = x(k) / A(k,k);
end
x = x';