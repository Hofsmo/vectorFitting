function [realPoles, complexPoles] = findPoles(x, y, t, complexPoles, realPoles, tol)
% FINDPOLES find poles of a system using vector fitting
%
% INPUT:
%   x: Input signal
%   y: Output signal
%   t: time signal
%   complexPoles: The initial complex poles of the system. The code
%   calculates the conjugates
%   realPoles: The initial real poles of the system
%
% OUTPUT:
%   poles: The estimated poles

if nargin < 6
    tol = 1e-5;
end

if nargin < 5
    realPoles = [];
end

% Number of complexPoles
nC = numel(complexPoles);

% Number of real Poles
nR = numel(realPoles);

% Timesteps
ts = numel(t);

% Initialize some variables for convenience
kn = 0; % Residues for realPoles
knI = 0; % Real part of the residues for complexPoles
knII = 0; % Imaginary part of the residues for complexPoles

xnR = sparse(ts,nR);
ynR = sparse(ts,nR);

xnC = sparse(ts,nC);
ynC = sparse(ts,nC);

% Convolution between exponential of each pole and signals. Results are
% stored in separate columns for each pole
if ~isempty(realPoles)
    xnR = windowConv (x, realPoles, t, false);
    ynR = windowConv (y, realPoles, t, false);
end

if ~isempty(complexPoles)
    xnC = windowConv (x, complexPoles, t, false);
    ynC = windowConv (y, complexPoles, t, false);  
end

A = [x, xnR, -ynR, real(xnC), -imag(xnC), -real(ynC), imag(ynC);...
     sparse(ts,1), sparse(ts, nR), sparse(ts,nR), imag(xnC), real(xnC), -imag(ynC), -real(ynC)];

sol = A\[y;sparse(ts,1)];

kn = sol(nR+2:1+2*nR)'; % Residues for real poles
knI = sol(end-2*nC+1:end-nC)'; % Real part of resideus for complex poles
knII = sol(end-nC+1:end)'; % Imagiary part of residues for complex poles

% Check whether or not we already have the correct poles
% Remember to check this assumption properly
if all(kn<tol) && all(knI<tol)
    return
%Check if we are dealing with complex pairs
elseif nC > 0
%Create the Â matrix from Gustavsen paper
    AHat = diag(real(complexPoles));
    %Indices of the superdiagonal
    v = 1:2:n-1;
    w = 2:2:n;
    superDiag = sub2ind([n,n], v, w); % Not really a superdiagonal
    AHat(superDiag) = abs(a);
    subDiag = sub2ind([n,n], w, v);
    AHat(subDiag) = -a;
    bc = zeros(n);
    bc(superDiag) = 2*imag(kn(1:2:end));
    idx = sub2ind([n,n], 1:2:n-1,1:2:n-1);
    bc(idx) = 2*real(kn(1:2:end));
    poles = eig(AHat-bc)';
elseif nR > 0
    %Find the zeros of the fit function
        realPoles=eig(diag(realPoles)-ones(numel(kn),1)*kn)';
else
    error('No poles specified');
end
% Flip unstable poles
%poles(real(poles)>0) = conj(poles(real(poles)>0)*-1);