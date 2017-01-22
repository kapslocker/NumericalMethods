A = [9,6,8; 6, 8 ,6; 8 , 6, 14];
%A = [9,1,2;1,1,5;2,5,7];
A = A*A';
b = [1;2;3];

%Cholesky_decomposition(A);
%L = chol(A,'lower')
X = [0:0.01:1];
Y = sin(2*3.14159*X);
%Y = sin(X)./X;
U = [0:0.005:1];
%Newton_DD(X,Y,U)
%V = piecewise_linear(X,Y,U);
V = cubic_spline(X,Y,U);
%Gauss(A,b)
%CG(A,b,A,[0;0;0])
%x = SOR(A,b,[0;0;0])
%A\b