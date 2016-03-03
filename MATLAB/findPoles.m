function [realPoles, complexPoles]...
    = findPoles(x, y, t, complexPoles, realPoles, tol)
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
%   realPoles: The estimated real poles
%   complexPoles: The estimated complex poles

if nargin < 6
    tol = 1e-4;
end

if nargin < 5
    realPoles = [];
end

% Timesteps
ts = numel(t);

% Initialize some variables for convenience
kn = 0; % Residues for realPoles
knI = 0; % Real part of the residues for complexPoles
knII = 0; % Imaginary part of the residues for complexPoles

poleForward = [];

def = true;

while def

    % Number of complexPoles
    nC = numel(complexPoles);

    % Number of real Poles
    nR = numel(realPoles);
    
    xnR = sparse(ts,nR);
    ynR = sparse(ts,nR);

    xnI = sparse(ts,nC);
    ynI = sparse(ts,nC);
    xnII = sparse(ts,nC);
    ynII = sparse(ts,nC);


    % Convolution between exponential of each pole and signals. Results are
    % stored in separate columns for each pole
    if ~isempty(realPoles)
        xnR = windowConv (x, realPoles, t);
        ynR = windowConv (y, realPoles, t);
    end

    if ~isempty(complexPoles)
        temp = windowConv (x, complexPoles, t);
        xnI = real(temp);
        xnII = imag(temp);

        temp = windowConv (y, complexPoles, t);
        ynI = real(temp);
        ynII = imag(temp);
    end

    A = [x, xnR, -ynR, 2*xnI, -2*xnII, -2*ynI, 2*ynII];

    sol = (full(A)\y);
    
    kn = sol(nR+2:1+2*nR)'; % Residues for real poles
    knI = sol(end-2*nC+1:end-nC)'; % Real part of residues for complex poles
    knII = sol(end-nC+1:end)'; % Imagiary part of residues for complex poles
    
    % Find the highest residue
    relkn = [abs(kn),abs(complex(knI,knII))];
    [~, iMaxkn] = max(relkn);
    
    % If no residues are zero check for large relative differences in
    % residue size
    minR = min (relkn);
    if minR
        relkn = min(relkn)/max(relkn);
    else
        relkn = 1;
    end
    
    % Check if the system is rank deficient or relative difference smaller
    % than tol
    if rank(full(A)) < 2*nR+4*nC+1 || relkn < tol
        def = true;
        if iMaxkn > nR
            complexPoles = complexPoles(knI~=knI(iMaxkn-nR));            
        else
            realPoles = realPoles(kn~=kn(iMaxkn));
        end
    else
        def = false;
    end
end
    
% Check whether or not we already have the correct poles
if all(abs(kn)<tol) && all(abs(knI) <tol)
    return
end
%Check if we are dealing with complex pairs
if nC > 0
%Create the Â matrix from Gustavsen paper
    temp = repmat(real(complexPoles.'),1,2)';
    AHat = diag(temp(:));
    %Indices of the superdiagonal
    v = 1:2:2*nC-1;
    w = 2:2:2*nC;
    superDiag = sub2ind([2*nC,2*nC], v, w); % Not really a superdiagonal
    AHat(superDiag) = imag(complexPoles);
    subDiag = sub2ind([2*nC,2*nC], w, v);
    AHat(subDiag) = -imag(complexPoles);
    bc = zeros(2*nC);
    bc(superDiag) = 2*knII;
    idx = sub2ind([2*nC,2*nC], 1:2:2*nC-1,1:2:2*nC-1);
    bc(idx) = 2*knI;
    temp = eig(AHat-bc)';
    
    % If there are real poles forward them
    poleForward = -abs(real(temp(abs(imag(temp))<tol)));
    temp = temp(abs(imag(temp))>tol);
    
    % Return only the negative pair of the complex conjugate pairs
    complexPoles = complex(-abs(real(temp(1:2:end))), -abs(imag(temp(1:2:end))));
end
if nR > 0
    % Find the zeros of the fit function
     realPoles=eig(diag(realPoles)-ones(numel(kn),1)*kn)';
    % See if there are any complex poles
    if ~isreal(realPoles)
        % Order the complex poles into pairs
        temp = cplxpair(realPoles(abs(imag(realPoles))>tol));
        % Put the poles not already in complexPoles into complexPoles
        complexPoles = [complexPoles,...
            setdiff(complex(-abs(real(temp(1:2:end))),...
        -abs(imag(temp(1:2:end)))),complexPoles)];
        realPoles = -abs(realPoles(abs(imag(realPoles))<tol));
    else
        % Flip unstable poles
        realPoles(realPoles>0) = -realPoles(realPoles>0);
    end
end
temp = setdiff(poleForward,realPoles);
if ~isempty(temp)
    % In case the code for complex poles found any real poles.
    realPoles = [realPoles, temp];
end
