function [zerr, perr] = realPolesTest ()
% Create the time vector
ts = 1e-3;
Tf = 5;
t = [0:ts:Tf];

% Create a pulse as input signal
x = (exp(-t)-exp(-2*t))'; 

% Create a superposition of decaying exponentials giving real poles as
% output
y = 1/2*(exp(-t)-exp(-3*t))';

% Create the initial poles
initPoles = -linspace(0,3,5);

[pn,cn,d] = fitVectorTime(x, y, t, initPoles,true);

zAnal = [-2]; % Analytical answer if x was a real step
pAnal = [-3]; % Analytical answer if x was a real step

zn = roots(residue(cn(abs(cn)>1e-3),pn(abs(cn)>1e-3),d));

zerr = zAnal-zn;
perr = pAnal-pn(abs(cn)>1e-3);