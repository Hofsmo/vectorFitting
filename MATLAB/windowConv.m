function wave = windowConv (x, poles, t, doTrapz)
% WINDOWCONV calculates the convolution from the paper

if nargin < 4
    doTrapz = false;
end

% Number of timesteps
ts = numel(t);

% Time window
dt = t(2)-t(1);
n = numel (poles);

if ~doTrapz
    temp = reshape(cell2mat(arrayfun(@(pole) conv(exp(pole*t), x), poles,...
        'UniformOutput', false)),[],n);
    wave = temp(1:ts,:)*dt;
else
    wave = reshape(cell2mat(arrayfun(@(pole) trapzConv(x,pole,t), poles,...
        'UniformOutput', false)),[],n);
end
