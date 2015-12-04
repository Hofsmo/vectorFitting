function poles = findPoles(x, y, t, initPoles, fReal)
% FINDPOLES find poles of a system using vector fitting
%
% INPUT:
%   x: Input signal
%   y: Output signal
%   t: time signal
%   initPoles: The initial poles of the system
%
% OUTPUT:
%   poles: The estimated poles

if nargin < 5
    fReal = false;
end
% Number of poles
n = numel(initPoles);

% Sort the poles
initPoles = sort(initPoles);

% Convolution between exponential of each pole and signals. Results are
% stored in separate columns for each pole
xn = windowConv (x, initPoles, t, false);
yn = windowConv (y, initPoles, t, false);

A = [x, xn, -yn];

sol = A\y;

kn = sol(end-n+1:end)';

poles = initPoles;

%Check if we are dealing with complex pairs
if any(imag(initPoles)>0 & ~fReal)
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



