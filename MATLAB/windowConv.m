function wave = windowConv (x, poles, t)
% WINDOWCONV calculates the convolution from the paper
%
% INPUT:
%   x: input signal to be convolved with the poles
%   poles: the poles to use in the convolution
%   t: the time vector
%
% OUTPUT:
%    wave: the resulting waveform from the convolution

% Number of timesteps
ts = numel(t);

% Time window
dt = t(2)-t(1);
n = numel (poles);

temp = reshape(cell2mat(arrayfun(@(pole) conv(exp(pole*t), x), poles,...
    'UniformOutput', false)),[],n);
wave = temp(1:ts,:)*dt;

