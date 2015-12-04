function [zerr,perr] = complexPolesTest ()

    ts = 0.001;
    Tf = 10;
    t = 0:ts:Tf;
    
    i = (50*exp(-0.6*t).*sin(0.8*t))';
    
    alpha = 100000;
    
    x = 160*(1-exp(-alpha*t))';
    
    betha = linspace(0.1,5,2);
    p = complex(-betha/100,betha);
    initPoles = [p,conj(p)];
    
    [zn, pn] = fitVectorTime(x, i, t, initPoles);
  
    step(tf(real(zn),pn),Tf);
   
    z = [0.25, 0]; % Analytical answer if x was a real step
    p = [1, 1.2, 1]; % Analytical answer if x was a real step
    zerr = z-zn;
    perr = p-pn;