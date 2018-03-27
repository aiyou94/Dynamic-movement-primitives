function [ts,xs]=canonicalsystem(steps,dt);

%constants
dt=0.001;
tsim=1;
steps=tsim/dt;
x0=1;
xf=0;
ax=8;

%euler integration to compute the timing signal
for n=1:steps;
    
    ts(1)=0;
    xs(1,1)=x0;
    xs(2,1)=-ax;
    ts(n+1)=ts(n)+dt;
    xs(1,n+1)=xs(1,n)+dt*xs(2,n);
    xs(2,n+1)=-ax*xs(1,n);
end

end

