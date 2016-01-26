function [err] = complexPolesTest ()

    ts = 0.001;
    Tf = 10;
    t = 0:ts:Tf;
    
    H = tf(conv([1,1],[1,2]),conv([1,1/2+50i],[1,1/2-50i]));
    x = (exp(-t)-exp(-2*t))'; 
    y = lsim(H,x,t);
    
    [p,z]=pzmap(H);
    
    betha = -linspace(45,62,3);
    initPoles = complex(betha/100,betha);
    
    [pn,cn,d] = fitVectorTime(x, y, t, initPoles);
    
    [den,num] = residue(cn(abs(cn)>1e-3),pn(abs(cn)>1e-3),d);
    yr = lsim(tf(den,num),x,t);
    err = immse(y,yr);