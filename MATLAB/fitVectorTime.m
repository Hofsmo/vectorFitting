function [pn,cn,d] = fitVectorTime(x, y, t, complexPoles, realPoles, directCoupling,...
    tol, i_max)
% FITVECTORTIME finds transfer function of a system using vector fitting
%
% INPUT:
%   x: Input signal
%   y: Output signal
%   t: time signal
%   initPoles: The initial poles of the system
%   directCoupling: Boolean to turn on or off the direct coupling term. The
%   default is false
%   fReal: Forces the part of the code not dealing with complex conjugate
%   pairs to be run.
%   tol: If the residues change less than this, we assume convergence
%   i_max: The maximum number of iterations
%
% OUTPUT:
%   cn: residues of the system
%   pn: poles of the system
%   d: direct coupling of the system

if nargin < 8
    i_max = 100;
end

if nargin < 7
    tol = 1e-5;
end

if nargin < 6
    directCoupling = true;
end

if nargin < 5
    realPoles = [];
end

error = 10;
i = 0;

while error > tol && i < i_max
    poles = findPoles(x, y, t, complexPoles, realPoles, tol);
    error = norm (poles-initPoles);
    initPoles = poles;
    i = i+1;
end

xn = windowConv (x, poles, t, false);

if directCoupling
    A = [x, xn];
    H = A\y;
    d = H(1);
    cn = H(2:end);
else
    cn = xn\y;
    d = 0;
end
pn = poles;
