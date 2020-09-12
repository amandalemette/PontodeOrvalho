clear all
clc
global P bnz tol y

y(1) = 0.5; % bezeno
y(2) = 0.5; % tolueno

P = 760; %mmHg

p = Antoine.data;
disp(p)

bnz = p.bnz;
tol = p.tol;

T0 = 90; % C
Torv = fzero(@fun, T0);

disp(Torv)

% Ki = Pisat/P = yi/xi --> xi = yi*P/Pisat

x(1) = y(1)*P/bnz.Psat(Torv);
x(2) = y(2)*P/tol.Psat(Torv);

disp(x)





function err = fun(T)
    global P bnz tol y
    err = 1/P - y(1)/bnz.Psat(T) - y(2)/tol.Psat(T);
end



