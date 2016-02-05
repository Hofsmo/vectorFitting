function wave = windowConv (x, poles, t)
% WINDOWCONV calculates the convolution from the paper

% Number of timesteps
ts = numel(t);

% Time window
dt = t(2)-t(1);
n = numel (poles);

temp = reshape(cell2mat(arrayfun(@(pole) conv(exp(pole*t), x), poles,...
    'UniformOutput', false)),[],n);
wave = temp(1:ts,:)*dt;

