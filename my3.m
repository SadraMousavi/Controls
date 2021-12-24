clc;clear all;

%%  Discrete MIMO SERVO + Observer
G=[1 -0.5 0.5;0 1 -1;1 0 0.7]
H=[1 0;0 0.5;0.5 0]
C=[1 0 0;0 1 0] % y1=x1 y2=x2
M=[H G*H]
rank(M)
N=[C;C*G]
rank(N)
Gh=[G zeros(3,2);-C eye(2)]
Hh=[H;zeros(2)]
Mh=[Hh Gh*Hh Gh^2*Hh Gh^3*Hh]
rank(Mh)
mu_c=[0.5+0.5*j 0.5-0.5*j -0.5 0.5 0.3];
Kc=place(Gh,Hh,mu_c)
mu_o=[0.02 -0.02 0.025]
Ke=place(G',C',mu_o)'

yr1 = @(t) 5*sign(sin(0.5*t));
yr2 = @(t) 10*sign(sin(0.25*t));
inp.K  = Kc;
inp.Ko = Ke;
inp.A  = G;
inp.B  = H;
inp.C  = C;
inp.y_ref1 = yr1;
inp.y_ref2 = yr2;
k=1;
h=0.05;
N=15;
X0=[1 -2 3 0 0]';
X0e=[0 0 0]';
Time = 1:h:N;

[X, Xh] = myRungeKutta2D(@Fun3, Time, X0, X0e, inp);
figure;
subplot(3,1,1);plot(Time,X(1,:),'b',Time,Xh(1,:),'r',Time,yr1(Time),'g');
subplot(3,1,2);plot(Time,X(2,:),'b',Time,Xh(2,:),'r',Time,yr2(Time),'g');
subplot(3,1,3);plot(Time,X(3,:),'b',Time,Xh(3,:),'r');


%% Continuous time with Digital Control
A=[1 2;-2 -3];
B=[1;-2];
C=[1 0];
h=0.3;
mu1c=-3;mu2c=-2;
mu1d=exp(mu1c*h);mu2d=exp(mu2c*h);
G=expm(A*h);
H=inv(A)*(expm(A*h)-eye(2))*B;
%%%% H=h*(eye(2)+A*h/2+A^2*h^2/6+...)*B
Kd=acker(G,H,[mu1d,mu2d]);% Kd=acker(G,H,0*[mu1d,mu2d]);
Kc=place(A,B,[mu1c,mu2c]);

T=10;
dt=0.01;
X0=[5;-10];
Time = 1:dt:T;


inp.K  = Kd;
% inp.K  = Kc;
inp.A  = A;
inp.B  = B;
inp.C  = C;
inp.h  = h;
inp.dt = dt;



[u,X] = myRungeKutta1N(@Fun5, Time, X0, inp);
figure;
subplot(3,1,1);plot(Time,X(1,:));
subplot(3,1,2);plot(Time,X(2,:));
subplot(3,1,3);plot(Time(1:end-1),u);




    
%% Functions

function [u,X] = myRungeKutta1N(Fun, tspan, X0, inp)
    dt = tspan(2)-tspan(1);
    X(:,1)  = X0;
    k=0;
    for i=1:length(tspan)-1
        t  = tspan(i);
        Xi = X(:,i);
        if mod(i,floor(inp.h/inp.dt))==1
            k=k+1;
            ud(k)=-inp.K*Xi;
        end
        u(i)=ud(k);
        K1 = Fun(t,Xi,inp,u(i));
        K2 = Fun(t+dt/2,Xi+K1*dt/2,inp,u(i));
        K3 = Fun(t+dt/2,Xi+K2*dt/2,inp,u(i));
        K4 = Fun(t+dt,Xi+K3*dt,inp,u(i));
        Xi = Xi+(K1+2*K2+2*K3+K4)/6*dt;
        X(:,i+1)=Xi;
    end
end

function [X, Xh] = myRungeKutta2D(Fun, tspan, X0, Xh0, inp)
    X(:,1)  = X0;
    Xh(:,1) = Xh0;
    for i=1:length(tspan)-1
        t  = tspan(i);
        Xi = X(:,i);
        Xhi = Xh(:,i);
        [X(:,i+1), Xh(:,i+1)] = Fun(t,Xi,Xhi,inp);
    end
end



function [dX,dXh]=Fun3(t,X,Xh,inp)
    yr1 = inp.y_ref1(t);
    yr2 = inp.y_ref2(t);
    y   = inp.C*X(1:end-2);
    u   = -inp.K*[Xh;X(end-1);X(end)];
    dX  = [inp.A*X(1:end-2)+inp.B*u ; yr1-y(1)+X(end-1); yr2-y(2)+X(end)];%moadele servo neveshte shavad
    
    yh  = inp.C*Xh;
    dXh = inp.A*Xh+inp.B*u+inp.Ko*(y-yh);
end


function dX=Fun5(t,X,inp,u)
    dX = inp.A*X+inp.B*u;
end