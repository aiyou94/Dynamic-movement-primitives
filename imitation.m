%input center & variance of kernels, & choice of weight vector
%w=[0 0 0 0 0 0 0 0 0 0];
%w=[0.1 0.2 0.3 0.4 .5 .6 .7 .8 .9 1];
%w=[1 5 10 20 50 20 10 5 2 1];
%w=[500 100 50 25 10 5 2 2 1 1];
%w=[1 2 3 4 5 6 7 8 9 10];
%w=[10000 100 1 1 1 1 1 1 1000 1000];
w=[10000 100 1 1 1 1 1 1 1000 1000];
center=[1 0.6294 0.3692 0.2494 0.1569 0.0988 0.0622 0.0391 0.0246 0.0155];
sig=[41.6667 16.3934 6.5459 2.5840 1.0235 0.4054 0.1606 0.0636 0.0252 0.0252]/1000;

%call the functions to get time(ts), phase(xs), & basis functions(psi) data
[ts,xs]=canonicalsystem(1000,0.001);
psi=basis(center,sig,xs);
[ys,yds,ydds]=dmp(w,psi);

%constants
T=1001;
g=1;
y0=0;
dt=0.001;
az=25;
bz=25;

%load data
load imitation.data;
y_demo=imitation(:,1);
yd_demo=imitation(:,2);
ydd_demo=imitation(:,3);

%goal and start of training data
goal=y_demo(1001,1);
start=y_demo(1,1);

% Compute the phi matrix
for idim=1:10
    for tt=1:T-1
    topp(tt+1,idim)=psi(tt,idim).*xs(1,tt).*(goal-start);
    bottomm=sum(psi,1);
end
end
phii=topp./bottomm;

%approximator of forcing function
f_target=ydd_demo-az.*(bz.*(goal-y_demo)-yd_demo);

%Using batch regression, compute the weight vector for the training data
w_learned=inv(phii'*phii)*((phii')*f_target);   %10x1
w_learned=w_learned';                        %1x10

%use the weight vector to compute the desired imitation trajectory position, velocity and acceleration

for idim=1:10
    for tt=1:T-1
    %for part c plots
    weighted=sum(psi.*repmat(w(:,:),T,1),2); 
    bottom=sum(psi,2);
    phi=(weighted./bottom);
    %
    weighted_learned=sum(psi.*repmat(w_learned(:,:),T,1),2);
    bottom_learned=sum(psi,2);
    phi_learned=(weighted_learned./bottom_learned);
    %epsilon(tt+1)=(xs(1,:)').*(goal-start);
    f_learned=phi_learned.*(xs(1,:)').*(goal-y_demo);
    ys_learned(1)=start;
    yds_learned(1)=-2;    %-2 was chosen based on visual inspection of first iteration
    % euler integration (two approaches to approximate forcing fxn)
    ydds_learned(tt+1)=az*(bz.*(g-ys_learned(tt))-yds_learned(tt))+f_target(tt);
    %ydds_learned(tt+1)=az*(bz*(goal-ys(tt))-yds(tt))+f_learned(tt);
    yds_learned(tt+1) = yds_learned(tt) + dt*ydds_learned(tt+1);
    ys_learned(tt+1)  =  ys_learned(tt) + dt* yds_learned(tt+1);
    % Maybe smooth by applying first order filter on g
    %g(1,1)=1;
    %g(2,1)=[az/8];
    %g(2,tt+1)=[az/8]*[g(1,tt)];
    %g(1,tt+1)=g(1,tt)+dt*g(2,tt);
    end    
end


%plots
figure(1)
subplot(2,1,1)
plot(ts,xs(1,:))
ylabel('x')
title('phase variable vs time')
subplot(2,1,2)
plot(ts,xs(2,:))
title('xd vs time')
xlabel('Time(s)')
ylabel('xd')

figure(2)
subplot(4,1,1)
stem(w)
title('weights vector')
subplot(4,1,2)
plot(ts,psi)
xlabel('Time(s)')
ylabel('psi(i)')
title('basis functions vs time')
subplot(4,1,3)
plot(xs(1,:),phi)
xlabel('x')
ylabel('psi(i)')
title('weighted & normalized summation vs phase variable(x)')
subplot(4,1,4)
plot(ts,phi)
title('weighted & normalized summation vs time(s)')
xlabel('Time(s)')
ylabel('phi(i)')

figure(3)
subplot(3,1,1)
plot(ts,ys,'k')
ylabel('y')
title('Position vs time')
subplot(3,1,2)
plot(ts,yds,'k')
title('Velocity vs time')
ylabel('yd')
subplot(3,1,3)
plot(ts,ydds,'k')
title('Acceleration vs time')
xlabel('Time(s)')
ylabel('ydd')


figure(4)
subplot(3,1,1)
plot(ts,ys_learned,'k',ts,y_demo,'r')
ylabel('y')
title('Position vs time')
subplot(3,1,2)
plot(ts,yds_learned,'k',ts,yd_demo,'r')
title('Velocity vs time')
ylabel('yd')
subplot(3,1,3)
plot(ts,ydds_learned,'k',ts,ydd_demo,'r')
title('Acceleration vs time')
xlabel('Time(s)')
ylabel('ydd')
legend('Learned data','training data');
