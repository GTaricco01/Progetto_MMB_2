clearvars;close all;clc
v = my_trap(@(t)t.^2,-1,1,100)



function v = my_trap(f,a,b,n)
% Composite Trapezoidal Rule
% a starting point
% b ending point
% n number of subintervals
v = (b-a)/n*(f(a)/2+sum(f(a+(1:(n-1))*(b-a)/n))+f(b)/2);
end
