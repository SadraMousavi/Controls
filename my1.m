clc;clear all;close all;

%% Regulator
% DX=AX+Bu,  y=CX

A=[1 2 -1;0 1 2;2 1 0]
B=[1 2 -1]'
C=[1 0 0]
M=[B A*B A^2*B]
rank(M)
N=[C;C*A;C*A^2]
rank(N)
eig(A)
mud=[-5+j -5-j -10]
pA=poly(A)
W=[pA(3) pA(2) 1;pA(2) 1 0;1 0 0]
T=M*W
pcl=conv(conv([1 -mud(1)],[1 -mud(2)]),[1 -mud(3)])
K=[pcl(4)-pA(4) pcl(3)-pA(3) pcl(2)-pA(2)]*inv(T)
eig(A-B*K)
K=place(A,B,mud)
K=[0 0 1]*inv(M)*(pcl(1)*A^3+pcl(2)*A^2+pcl(3)*A+pcl(4)*eye(3))
K=acker(A,B,mud)

inp.K = K;
inp.A = A;
inp.B = B;
inp.C = C;
% inp.D = D;
T=10;
dt=0.05;
Time = 0:dt:T;
X0=[3;-3;5];
X = myRungeKutta1(@Fun1, Time, X0, inp);
figure;
subplot(3,1,1);plot(Time,X(1,:));
subplot(3,1,2);plot(Time,X(2,:));
subplot(3,1,3);plot(Time,X(3,:));


%% SERVO

Ah=[A zeros(3,1);-C 0]
Bh=[B;0]
Mh=[Bh Ah*Bh Ah^2*Bh Ah^3*Bh]
rank(Mh)
mud=[-10 -15 -17 -20]
K=place(Ah,Bh,mud)


yr = @(t) 5*sign(sin(0.5*t));
inp.K = K;
inp.A = A;
inp.B = B;
inp.C = C;
inp.y_ref = yr;
T=20;
dt=0.01;
Time = 0:dt:T;
X0=[0;0;0;0];
X = myRungeKutta1(@Fun2, Time, X0, inp);
figure;
subplot(3,1,1);plot(Time,X(1,:),'b',Time,yr(Time),'r');
subplot(3,1,2);plot(Time,X(2,:));
subplot(3,1,3);plot(Time,X(3,:));


%% SERVO+Full Order Observer
%%% DXh=A*Xh+B*u+Ko*(y-yh)

N=[C;C*A;C*A^2]
rank(N)
muo=[-10 -10 -10];
Ko=acker(A',C',muo);
Ko=Ko';
eig(A-Ko*C);
Ah=[A zeros(3,1);-C 0]
Bh=[B;0]
Mh=[Bh Ah*Bh Ah^2*Bh Ah^3*Bh]
rank(Mh)
mud=[-2+j -2-j -3 -3]
K=acker(Ah,Bh,mud)

yr = @(t) 5*sign(sin(0.5*t));
inp.K  = K;
inp.Ko = Ko;
inp.A  = A;
inp.B  = B;
inp.C  = C;
inp.y_ref = yr;
T=20;
dt=0.01;
Time = 0:dt:T;
X0=[0;0;0;0];
Xh0=[1;1;-1];

[X, Xh] = myRungeKutta2(@Fun3, Time, X0, Xh0, inp);
figure;
subplot(3,1,1);plot(Time,X(1,:),'b',Time,yr(Time),'r',Time,Xh(1,:),'g');
subplot(3,1,2);plot(Time,X(2,:),Time,Xh(2,:),'g');
subplot(3,1,3);plot(Time,X(3,:),Time,Xh(3,:),'g');



%% SERVO (Feed Forward) +Full Order Observer
%%% DXh=A*Xh+B*u+Ko*(y-yh)

N=[C;C*A;C*A^2]
rank(N)
muo=[-10 -10 -10];
Ko=acker(A',C',muo);
Ko=Ko';
eig(A-Ko*C);
mud=[-2+j -2-j -3]
K=acker(A,B,mud)

yr = @(t) 5*sign(sin(0.5*t));
inp.K  = K;
inp.Ko = Ko;
inp.A  = A;
inp.B  = B;
inp.C  = C;
inp.y_ref = yr;
T=20;
dt=0.01;
Time = 0:dt:T;
X0=[0;0;0];
Xh0=[1;1;-1];

[X, Xh] = myRungeKutta2(@Fun4, Time, X0, Xh0, inp);
figure;
subplot(3,1,1);plot(Time,X(1,:),'b',Time,yr(Time),'r',Time,Xh(1,:),'g');
subplot(3,1,2);plot(Time,X(2,:),Time,Xh(2,:),'g');
subplot(3,1,3);plot(Time,X(3,:),Time,Xh(3,:),'g');



%% Functions
function [Time,X] = myRungeKutta(Fun, tspan, X0)
    T  = tspan(end);
    dt = tspan(2)-tspan(1);
    t  = tspan(1);
    X(:,1)  = X0;
    Time(1) = t;
    i = 1;
    while t < T
        Xi=X(:,i);
        K1=Fun1(t,Xi,u);
        K2=Fun1(t+dt/2,Xi+K1*dt/2,u);
        K3=Fun1(t+dt/2,Xi+K2*dt/2,u);
        K4=Fun1(t+dt,Xi+K3*dt,u);
        Xi=Xi+(K1+2*K2+2*K3+K4)/6*dt;
        X(:,i+1)=Xi;
        Time(i+1)=t+dt;
        i=i+1;
        t=t+dt;
    end
end

function X = myRungeKutta1(Fun, tspan, X0, inp)
    dt = tspan(2)-tspan(1);
    X(:,1)  = X0;
    for i=1:length(tspan)-1
        t  = tspan(i);
        Xi = X(:,i);
        K1 = Fun(t,Xi,inp);
        K2 = Fun(t+dt/2,Xi+K1*dt/2,inp);
        K3 = Fun(t+dt/2,Xi+K2*dt/2,inp);
        K4 = Fun(t+dt,Xi+K3*dt,inp);
        Xi = Xi+(K1+2*K2+2*K3+K4)/6*dt;
        X(:,i+1)=Xi;
    end
end

function [X, Xh] = myRungeKutta2(Fun, tspan, X0, Xh0, inp)
    dt = tspan(2)-tspan(1);
    X(:,1)  = X0;
    Xh(:,1) = Xh0;
    for i=1:length(tspan)-1
        t  = tspan(i);
        Xi = X(:,i);
        Xhi = Xh(:,i);
        [D1, K1] = Fun(t,Xi,Xhi,inp);
        [D2, K2] = Fun(t+dt/2,Xi+D1*dt/2,Xhi+K1*dt/2,inp);
        [D3, K3] = Fun(t+dt/2,Xi+D2*dt/2,Xhi+K2*dt/2,inp);
        [D4, K4] = Fun(t+dt,Xi+D3*dt,Xhi+K3*dt,inp);
        Xi = Xi+(D1+2*D2+2*D3+D4)/6*dt;
        X(:,i+1)=Xi;
        Xhi = Xhi+(K1+2*K2+2*K3+K4)/6*dt;
        Xh(:,i+1)=Xhi;
    end
end

function dX=Fun1(t,X,inp)
    u  = -inp.K*X;
    dX = inp.A*X+inp.B*u;
end

function dX=Fun2(t,X,inp)
    yr = inp.y_ref(t);
    y  = inp.C*X(1:end-1);
    u  = -inp.K*X;
    dX = [inp.A*X(1:end-1)+inp.B*u ; yr-y];
end

function [dX,dXh]=Fun3(t,X,Xh,inp)
    yr  = inp.y_ref(t);
    y   = inp.C*X(1:end-1);
    u   = -inp.K*X;
    dX  = [inp.A*X(1:end-1)+inp.B*u ; yr-y];
    
    yh  = inp.C*Xh;
    dXh = inp.A*Xh+inp.B*u+inp.Ko*(y-yh);
end

function [dX,dXh]=Fun4(t,X,Xh,inp)
    yr  = inp.y_ref(t);
    y   = inp.C*X;
    uff = -yr/(inp.C*inv(inp.A-inp.B*inp.K)*inp.B);
    u   = -inp.K*X+uff;
    dX  = [inp.A*X+inp.B*u];
    
    yh  = inp.C*Xh;
    dXh = inp.A*Xh+inp.B*u+inp.Ko*(y-yh);
end