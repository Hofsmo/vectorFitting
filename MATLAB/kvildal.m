
%data = xlsread('KvilldalT4_5normal.xlsx');

f = data(:,1);
f = f(~isnan(f));
u = data(:,2).*(cos(data(:,3)*pi/180) +1i*sin(data(:,3)*pi/180));
I = -data(:,4).*(cos(data(:,5)*pi/180) +1i*sin(data(:,5)*pi/180));
p = real (u.*conj(I))*3/10^6;
p = p(~isnan(p));
f = f(1:2200*50);
p = p(1:2200*50);


% wSize = numel(1:0.02:120);
% wts=ones(wSize,1)/wSize;
%  
% p1 =conv(p,wts,'valid');
% f1 = conv(f,wts,'valid');

% p1 = filter(Hlp,p);
% f1 = filter(Hlp,f);



%wSize = numel(1:0.02:120);
%   
% p1 = smooth(p,500);
% f1 = smooth(f,500);
t=[0:0.02:(length(f)-1)*0.02];
t = downsample(t,10);
%p1 = downsample(p1,10);
%f1= downsample(f1,10);

p1 = fnval(t,csaps(t,downsample(p,10),0.01));
f1 = fnval(t,csaps(t,downsample(f,10),0.01));

% f1=f1.';
% p1=p1.';
% 
f1 = mean(f1) -  f1;
p1 = p1 - mean(p1);

initPoles = -linspace(0.01,2*2*pi,20);
%initPoles = complex(initPoles/100,initPoles);
%initPoles = [initPoles,conj(initPoles)];

[den,num,pn,cn,d] = fitVectorTime (f1',p1',t, initPoles);