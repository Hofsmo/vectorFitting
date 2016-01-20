function y = trapzConv(x,p,t)

y = zeros(1,numel(t));
for i = 2:numel(t)
    y(i) = trapz(t(1:i), x(1:i).*exp(p*(t(i)-t(1:i)))');
end
