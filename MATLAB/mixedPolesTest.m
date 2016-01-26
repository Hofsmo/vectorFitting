function err = mixedPolesTest ()

T = 10;
ts = 1e-3;
t=0:ts:T;

x = (exp(-t)-exp(-2*t))';

H = tf(conv([1,2],[1,1]),conv(conv([1,1/2+50i],[1,1/2-50i]),[1,3]));

y = lsim(H,x,t);

realPoles = -linspace(0,5,3);
Betha = -linspace(42,60,3);
complexPoles = complex(Betha/100,Betha);

[pn,cn,d]=fitVectorTime(x,y,t,complexPoles,realPoles);

[den,num] = residue(cn(abs(cn)>1e-3),pn(abs(cn)>1e-3),d);
yr = lsim(tf(den,num),x,t);
err = immse(y,yr);