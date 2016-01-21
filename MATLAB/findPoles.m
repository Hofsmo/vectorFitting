function poles = findPoles(x, y, t, complexPoles, realPoles, tol)
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

% Convolution between exponential of each pole and signals. Results are
% stored in separate columns for each pole
if ~isempty(realPoles) && ~isempty(complexPoles)
    xnR = windowConv (x, realPoles, t, false);
    ynR = windowConv (y, realPoles, t, false);
    
    xnC = windowConv (x, complexPoles, t, false);
    ynC = windowConv (y, complexPoles, t, false);
    
    A = [x, xnR, -ynR, real(xnC), -imag(xnC), -real(ynC), imag(ynC);
         zeros(size(x)), zeros(size(xnR)), zeros(size(ynR)), imag(xnC), real(xnC), -imag(ynC), -real(ynC)];
elseif isempty(realPoles)
    xnC = windowConv (x, complexPoles, t, false);
    ynC = windowConv (y, complexPoles, t, false);
    
    A = [x, real(xnC), -imag(xnC), -real(ynC), imag(ynC);
         zeros(size(x)), imag(xnC), real(xnC), -imag(ynC), -real(ynC)];
     
elseif isempty(complexPoles)
    xnR = windowConv (x, realPoles, t, false);
    ynR = windowConv (y, realPoles, t, false);
    
     A = [x, xnR, -ynR];
else
    error('No poles specified')
end
% 



sol = A\[y; zeors(size(y))];
kn = sol(end-n+1:end)'; % Finding the residues

poles = initPoles;

% Check whether or not we already have the correct poles
% Remember to check this assumption properly
if all(real(kn)<tol)
    poles = initPoles;
%Check if we are dealing with complex pairs
elseif any(imag(initPoles)>0)
%Create the Â matrix from Gustavsen paper
    AHat = diag(real(initPoles));
    %Indices of the superdiagonal
    v = 1:2:n-1;
    w = 2:2:n;
    superDiag = sub2ind([n,n], v, w); % Not really a superdiagonal
    a = abs(imag(initPoles(1:2:end-1)));
    AHat(superDiag) = a;
    subDiag = sub2ind([n,n], w, v);
    AHat(subDiag) = -a;
    bc = zeros(n);
    bc(superDiag) = 2*imag(kn(1:2:end));
    idx = sub2ind([n,n], 1:2:n-1,1:2:n-1);
    bc(idx) = 2*real(kn(1:2:end));
    poles = eig(AHat-bc)';
else
    %Find the zeros of the fit function
        poles=eig(diag(initPoles)-ones(numel(kn),1)*kn)';
end
% Flip unstable poles
poles(real(poles)>0) = conj(poles(real(poles)>0)*-1);



