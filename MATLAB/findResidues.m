function [A,H] = findResidues(x,y,t, realPoles, complexPoles, directCoupling)
% FINDRESIDUES find the final residues
%
% Calculates the
% INPUT:
%
% x: input signal
% y: system respsonse
% t: time signal
% realPoles: real poles of the system
% complexPoles: complex conjugate poles of the system
% directCoupling: Boolean to turn on or off the direct coupling term. The
%   default is false
%
% OUTPUT.
% A: Signal matrix
% H: Residues

if nargin <6
    directCoupling = false;
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

