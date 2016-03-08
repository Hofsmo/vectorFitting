function wave = windowConv (x, poles, t)
% WINDOWCONV calculates the convolution from the paper
%
% Function that performs the convolution integral described in the original 
% TD-VF paper. However, using trapezoidal integration.
%
% INPUT:
%   x: The signal to convolve with the poles
%   poles: The poles to use in the convolution
%   t: the time signal
%
% OUTPUT:
%   wave: vector containing the waves resulting from the convolution
 
dt = t(2)-t(1); % Time step
n = numel (poles); % Number of poles 

wave = zeros(numel(t), n);

for i = 1:numel(n)
    % Constants for trapezoidal integration. They can be found in
    % paper comparing TD-VF and ZD-VF
    
	alpha = (1+poles(i)*dt/2)/(1-poles(i)*dt/2);
    lambda = dt/2/(1-poles(i)*dt/2);
    wave(:,i) = filter([lambda,lambda],[1,-alpha],x);
end


