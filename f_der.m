function val = f_der(x)
% Written by: Kapil Ahuja
%            (2014MT60663)
% This is required in the implementation of Newton's Method, to provide the
% gradient of the function defined in f.m .
% 
%%
val = 0.5*x.*sinh(x/4) + 2*cosh(x/4);