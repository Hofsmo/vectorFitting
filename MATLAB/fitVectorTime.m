function [pn,cn,d] = fitVectorTime(x, y, t, complexPoles, realPoles, directCoupling,...
    tol, i_max)
% FITVECTORTIME finds transfer function of a system using vector fitting
%
% INPUT:
%   x: Input signal
%   y: Output signal
%   t: time signal
%   complexPoles: The initial complex conjugate poles of the system. Only
%   provide the negative halves of the pairs. The poles should be given as
%   a row vector.
%   realPoles: The initial real poles of the system. They have to be
%   negative and given as a row vector.
%   directCoupling: Boolean to turn on or off the direct coupling term. The
%   default is false
%   tol: If the poles change less than this, we assume convergence
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
    tol = 1e-4;
end

if nargin < 6
    directCoupling = true;
end

if nargin < 5
    realPoles = [];
end

err = 10;
i = 0;

while err > tol && i < i_max
    [tempReal, tempComplex] =...
        findPoles(x, y, t, complexPoles, realPoles, tol);
    
    % Check whether or not the poles are moving
    if numel(tempReal)==numel(realPoles) &&...
            numel(tempComplex)==numel(complexPoles)
        err = norm(sort(tempReal) - sort(realPoles))...
            + norm(sort(tempComplex) - sort(complexPoles)); 
    end
            
    realPoles = tempReal;
    complexPoles = tempComplex;
    i = i+1;
end

if i ==i_max
    disp('iteration limit exceeded')
end

ts = numel(t); % Number of timesteps
nR = numel(realPoles);
nC = numel(complexPoles);

% Intialize variables to store the poles
xnR = sparse(ts, nR);
xnI = sparse(ts, nC);
xnII = sparse(ts, nC);

if ~isempty(realPoles)
    xnR = windowConv (x, realPoles, t);
end

if ~isempty(complexPoles)
    temp = windowConv (x, complexPoles, t);
    xnI = real(temp);
    xnII = imag(temp);
end

if ~directCoupling
    x = sparse(ts,1);
end
    
A = [2*xnI, -2*xnII, xnR, x];
     
H = full(A)\y; % It seems mldivide is not good enough for sparce matrices

d = full(H(end));
cnI = H(1:nC);
cnII = H(nC+1:2*nC);
cnR = H(2*nC+1:2*nC+nR);

% Check whether or not some of the residues are below the tolerance limit
idxR = cnR(abs(cnR)<tol);
idxC = cnI(abs(cnI)<tol) & cnII(abs(cnII)<tol);

if any(idxR) || any(idxC)
    [pn,cn,d] = fitVectorTime(x, y, t, complexPoles(~idxC),...
        realPoles(~idxR), directCoupling, tol, i_max);
else      
    pn = [realPoles, cplxpair([complexPoles, conj(complexPoles)])];
    cn = [full(cnR'), complex(full(cnI'), full(cnII')), conj(complex(full(cnI'), full(cnII')))];
end
