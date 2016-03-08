function [pn,cn,d] = fitVectorTime(x, y, t, complexPoles, realPoles, directCoupling,...
    lsq, tol, i_max)
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

if nargin < 9
    i_max = 100;
end

if nargin < 8
    tol = 1e-5;
end

if nargin < 7
    lsq = false;
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
    [tempReal, tempComplex, err] =...
        findPoles(x, y, t, complexPoles, realPoles, tol);
       
    if lsq
        [A,H]=findResidues(x,y,t, realPoles, complexPoles, directCoupling);
        temp=norm(y-A*H);
        if temp<err
            err = temp;
        end
    end
    realPoles = tempReal;
    complexPoles = tempComplex;
    i = i+1;
end

if i ==i_max
    disp('iteration limit exceeded')
end

[~,H] = findResidues(x,y,t, realPoles, complexPoles, directCoupling);

nR = numel(realPoles);
nC = numel(complexPoles);

d = full(H(end));
cnI = H(1:nC);
cnII = H(nC+1:2*nC);
cnR = H(2*nC+1:2*nC+nR);

if isempty(complexPoles)
    cnI = [];
    cnII = [];
    complexPoles = [];
end

if isempty(realPoles)
    cnR = [];
    realPoles = [];
end

% Check whether or not some of the residues are below the tolerance limit
idxR = abs(cnR')./abs(realPoles)<tol;
complexC = abs(complex(cnI,cnII))';
idxC = complexC./abs(complexPoles)<tol;

if any(idxR) || any(idxC)
    [pn,cn,d] = fitVectorTime(x, y, t, complexPoles(~idxC),...
        realPoles(~idxR), directCoupling, tol);
else      
    pn = [realPoles, cplxpair([complexPoles, conj(complexPoles)])];
    cn = [full(cnR'), complex(full(cnI'), full(cnII')), conj(complex(full(cnI'), full(cnII')))];
end