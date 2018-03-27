function [ys,yds,ydds]=dmp(w,psi);
%constants
az=25;
bz=25;    
g=1;
y0=0;
dt=0.001;
tsim=1;
steps=tsim/dt;
xs=canonicalsystem(steps,dt);
T=length(xs);

for idim=1:10
    for tt=1:T-1
    weighted=sum(psi.*repmat(w(:,:),T,1),2); 
    bottom=sum(psi,2);
    phi=(weighted./bottom);
    f=phi.*(xs(1,:)')*(g-y0);  %nonlinear fxn
    ys(1)=y0;
    yds(1)=0;
    %euler integration
    ydds(tt+1)=az*(bz*(g-ys(tt))-yds(tt))+f(tt);
    yds(tt+1) = yds(tt) + dt*ydds(tt+1);
    ys(tt+1)  =  ys(tt) + dt* yds(tt+1);
    end
end

end
