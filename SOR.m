function x = SOR(A,b,x0)
%Written by: Kapil Ahuja
%            (2014MT60663)
% Implements SOR with a parameter w = 1.5(arbitrarily chosen between 1-2).
% The iteration Matrix G is given by G = I - M_inv * A. (G = inv(M)*N, A = M-N,
% and the update conditions give this result.)
%%
s = size(A);
n = s(1);
E = zeros(n);
D = zeros(n);
for i = 1:n
    for j = 1:i
        E(i,j) = A(i,j);        %setting matrix E  = L + U.
    end
end
for i = 1:n
    D(i,i) = A(i,i);
end
w = 1.5;                        % The relaxation parameter, chosen arbitrarily.
x=x0;
r = b - A*x;
M_inv = w*inv(((1-w)*D + w*E)); % The matrix to be multiplied later in each iteration. 
i = 0;
while(1)
    r = b - A*x;                % evaluate error, for checking stopping criteria.
    if( i>100  || max(abs(r))<1e-6 ) % stop if absolute convergence, or if the solution diverges.
        m = sprintf('Solution found after %d iterations:',i);
        disp(m);
        break
    end
    x = x + M_inv*r;            % update solution.
    G = eye(n) - M_inv * A;
    if(max(eig(G))>1)           % Diverging condition is when spectral radius og G is greater than 1.
        disp 'Warning! The function will diverge for current value of w.'
        break
    end
    i = i+1;
end