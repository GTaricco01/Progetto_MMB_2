clearvars; close all;clc
N = 1e6;
f      = -10+9.1*rand(N,1);
beta   = .5;
gamma  = 0.024;

for i = 1:3
    f = MonteCarlo(f,beta,gamma,N,i);
end
