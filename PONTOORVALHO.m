clear all
clc

y(1) = 0.5; % benzeno
y(2) = 0.5; % tolueno

T = 113; % C

p = Antoine.data;
disp(p)

bnz = p.bnz;
tol = p.tol;

Pinv = y(1)/bnz.Psat(T) + y(2)/tol.Psat(T);
P = 1/Pinv;
disp(P)

% Ki = Pisat/P = yi/xi --> 1/xi = Pisat/(P yi) --> xi = (P yi)/Pisat

x(1) = P*y(1)/bnz.Psat(T);
x(2) = P*y(2)/tol.Psat(T);
disp(x)

