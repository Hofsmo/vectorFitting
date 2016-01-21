function [zerr,perr] = complexPolesTest ()

    ts = 0.001;
    Tf = 10;
    t = 0:ts:Tf;
    
    y = (105*exp(-t).*sin(1/105*t))';
    
    x = (exp(-t)-exp(-2*t))'; 
    
    betha = linspace(98,114,3);
    p = complex(-betha/100,betha);
    initPoles = [p,conj(p)];
    %p = roots([1,2,105^2])';
    
    [pn,cn,d] = fitVectorTime(x, y, t, initPoles);
   
    z = [-1, -2]; % Analytical answer if x was a real step
    p = roots([1,2,105^2])'; % Analytical answer if x was a real step
    
    zn = roots(residue(cn(abs(cn)>1e-3),pn(abs(cn)>1e-3),d))';
    zerr = z-zn;
    perr = p-pn;