function vf = runVecFit(x,y,t,complexPoles,realPoles,tol)
if nargin < 6
    tol =1e-5;
end
tic
[vf.pn,vf.cn,vf.d] = fitVectorTime (x, y, t, complexPoles, realPoles, false,...
    false,tol);
vf.time = toc;
[den,num] = residue(vf.cn(abs(vf.cn)>tol),vf.pn(abs(vf.cn)>tol),vf.d);
vf.fit =tf(real(den),num);
