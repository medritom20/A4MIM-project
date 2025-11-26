% This is just a testing script for development purposes.

clear
clc

A = load("bcsstk03.mat").Problem.A;

n = length(A);
m = 6;

x_ex = rand(n,m);
b = A * x_ex;
x_0 = zeros(n,m);

[x] = Algorithm4OLBCG(A,b,x_0,1e-8,200);

% semilogy(omeg) missing omeg

norm(x_ex - x)


A = load("s3dkt3m2.mat").Problem.A;

n = length(A);
m = 1;

x_ex = rand(n,m);
b = A * x_ex;
x_0 = zeros(n,m);
L = ichol(A,struct('type','ict','droptol',1e-5,'diagcomp',1e-2));

[x, omeg] = PBlockCG_OLeary(A,b,x_0,1e-8,3000,inv(L * L'));

figure(2)
semilogy(omeg)

norm(x_ex - x)