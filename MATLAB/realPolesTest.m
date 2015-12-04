function [zerr, perr] = realPolesTest ()
% Create the input function x
alpha = 1000; % The inverse of the time constant
ts = 1e-3;
Tf = 5;

t = [0:ts:Tf];
x = (1-exp(-alpha*t)'); % Approximation of a step function
y = t';

% Create the initial poles
initPoles = -linspace(0,3,5);

[zn,pn,cn,d] = fitVectorTime(x, y, t, initPoles,true);

[yfit, T] = step(tf(real(zn),pn),Tf);

plot(T,yfit);
hold
plot(T,T);

z = [0, 1]; % Analytical answer if x was a real step
p = [1, 0]; % Analytical answer if x was a real step
zerr = z-zn
perr = p-pn