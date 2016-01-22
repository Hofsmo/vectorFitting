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
    [tempReal, tempComplex] = findPoles(x, y, t, complexPoles, realPoles, tol);
    error = norm ([tempReal,tempComplex]-[realPoles,complexPoles]);
    realPoles = tempReal;
    complexPoles = tempComplex;
    i = i+1;
end

ts = numel(t); % Number of timesteps
nR = numel(realPoles);
nC = numel(complexPoles);

% Intialize variables to store the poles
xnR = sparse(ts, nR);
xnC = sparse(ts, nC);

if ~isempty(realPoles)
    xnR = windowConv (x, realPoles, t, false);
end

if ~isempty(complexPoles)
    xnC = windowConv (x, complexPoles, t, false);
end

if directCoupling
    A = [real(xnC), -imag(xnC), xnR, x;...
         imag(xnC), real(xnC), sparse(ts, nR), sparse(ts, 1)];
     
    H = A\[y;sparse(ts,1)];
    
    d = full(H(end));
    cnI = H(1:nC);
    cnII = H(nC+1:2*nC);
    cnR = H(2*nC+1:2*nC+nR);
else
    cn = xn\y;
    d = 0;
end
pn = full([realPoles, complexPoles, conj(complexPoles)]);
cn = full([cnR, complex(full(cnI), full(cnII)), conj(complex(full(cnI), full(cnII)))]);
