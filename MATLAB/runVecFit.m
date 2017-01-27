function vf = runVecFit(data,complexPoles,realPoles,tol)
if nargin < 4
    tol =1e-5;
end
tic
[vf.pn,vf.cn,vf.d] = fitVectorTime (data.InputData,...
    data.OutputData, data.SamplingInstants, complexPoles, realPoles, false,...
    false,tol);
vf.time = toc;
[den,num] = residue(vf.cn(abs(vf.cn)>tol),vf.pn(abs(vf.cn)>tol),vf.d);
vf.fit =tf(real(den),num);
