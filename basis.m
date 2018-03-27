function [psi]=basis(center,sig,xs);

%center and variances of the 10 basis functions
sig=[41.6667 16.3934 6.5459 2.5840 1.0235 0.4054 0.1606 0.0636 0.0252 0.0252]/1000;
center=[1 0.6294 0.3692 0.2494 0.1569 0.0988 0.0622 0.0391 0.0246 0.0155];

% number of basis functions (10 in this example)
n_basis=length(center);



%call canonical function
dt=0.001;
tsim=1;
steps=tsim/dt;
[ts,xs]=canonicalsystem(steps,dt);

%initialize array for psi
psi=zeros(length(xs),n_basis);

for i=1:n_basis;
     psi(:,i)=exp([-0.5/sig(i)].*(xs(1,:)-center(i)).^2);
end

end

