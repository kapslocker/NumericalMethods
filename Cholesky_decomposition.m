function L = Cholesky_decomposition(A)
%Written by: Kapil Ahuja
%            (2014MT60663)
% The condition for A to be symmetrix is checked first.
% At each step, if the generated diagonal element is not positive, then the
% matrix is not positive definite. The matrix is generated as the solution
% of the product: 
% [a11 a12 a13 a14 ...;     [g11                      [ g11 g21 g31 g41 .
%  a21 a22 a23 a24 ...;      g21 g22                        g22 g32 . . .
%  .                         g31 g32 g33            *           g33 . . .
%  .                    =                               
%  .                                                    
%  an1 an2 an3 an4 ...]      gn1 gn2 . . . . .gnn]                    gnn]
% 
% Each aii is then updated in place and the lower part is used as L. 
%%
s = size(A);
n = s(1);
B =A.';
if(~isequal(A,B))
    error('The Matrix is not symmetric.')
end
for k = 1:n-1
    if(A(k,k)<=0)
        error('The Matrix is not positive definite.')
    end
    A(k,k) = sqrt(A(k,k));
    for i = k+1:n
        A(i,k)=A(i,k)/A(k,k);
    end
    for j = k+1:n
        for i = j:n
            A(i,j) = A(i,j)-A(i,k)*A(j,k);
        end
    end
end
if(A(n,n)<=0)
    error('The matrix is not positive definite.')
end
%% Store only the lower part of A.
A(n,n) = sqrt(A(n,n));
for i=1:n
    for j=i+1:n
        A(i,j)=0;
    end
end
L = A;