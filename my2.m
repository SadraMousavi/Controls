clc;clear all;
%% EESA Regulator
A=[1 2 -1 0;0 2 3 1;-2 -1 0 1;1 1 -2 2]
B=[1 2;-1 0;0 1;2 -1]
C=[1 0 0 0;0 1 0 0]
M=[B A*B A^2*B]
rank(M)
lam1=-2;
lam2=-3;
lam3=-1+j;
lam4=-1-j;
S1=null([A-lam1*eye(4) B])
S2=null([A-lam2*eye(4) B])
S3=null([A-lam3*eye(4) B])
S4=null([A-lam4*eye(4) B])
%%%[v1;q1] [v2;q2] [v3;q3] [v4;q4]
v1q1=S1(:,1)
v2q2=S2(:,1)
v3q3=real(S3(:,1))
v4q4=imag(S3(:,1))
v1=v1q1(1:4);q1=v1q1(5:6);
v2=v2q2(1:4);q2=v2q2(5:6);
v3=v3q3(1:4);q3=v3q3(5:6);
v4=v4q4(1:4);q4=v4q4(5:6);
K=-[q1 q2 q3 q4]*inv([v1 v2 v3 v4])
eig(A-B*K)


inp.K  = K;
inp.A  = A;
inp.B  = B;
inp.C  = C;
T=20;
dt=0.01;
Time = 0:dt:T;
X0=[2;2;2;2];

[X] = myRungeKutta1(@Fun1, Time, X0, inp);
figure;
subplot(4,1,1);plot(Time,X(1,:));
subplot(4,1,2);plot(Time,X(2,:));
subplot(4,1,3);plot(Time,X(3,:));
subplot(4,1,4);plot(Time,X(4,:));


%% GCCF Regulator
Mh=[B(:,1) A*B(:,1) B(:,2) A*B(:,2)]  %gamma1=2, gamma2=2 gamma1>=gamma2
rank(Mh)
Mhi=inv(Mh)
e12=Mhi(2,:);
e22=Mhi(4,:);
Ti=[e12;e12*A;e22;e22*A]
T=inv(Ti)
% DX=AX+Bu---->DZ=AgZ+Bgv,    X=TZ, u=Fw,  w=v-HZ
% Ag=Ti*A*T-Ti*B*F*H   Bg=Ti*B*F
Ag=[0 1 0 0;0 0 0 0;0 0 0 1;0 0 0 0]
Bg=[0 0;1 0;0 0;0 1]
syms f1 f2 f3 f4
F=[f1 f2;f3 f4]
Bg-Ti*B*F  %  No solution for gamma1=3, gamma2=1
F=[1 0;0 1]
Bg-Ti*B*F
% Ti*A*T-Ag=Bg*H---->Bg'*(Ti*A*T-Ag)=H
H=Bg'*(Ti*A*T-Ag)
mu1=-1;mu2=-2;mu3=-2;mu4=-3;
% (s-mu1)*(s-mu4)=s^2+4*s+3
% (s-mu2)*(s-mu3)=s^2+4*s+4
Ad=[0 1 0 0;-3 -4 0 0;0 0 0 1;0 0 -4 -4]
% Ad=Ag-Bg*Gam---> Gam=Bg'*(Ag-Ad)
% K=-F*(H+Gam)*Ti
Gam=Bg'*(Ag-Ad)
K=F*(H+Gam)*Ti
eig(A-B*K)



inp.K  = K;
T=20;
dt=0.01;
Time = 0:dt:T;
X0=[2;2;2;2];

[X] = myRungeKutta1(@Fun1, Time, X0, inp);
figure;
subplot(4,1,1);plot(Time,X(1,:));
subplot(4,1,2);plot(Time,X(2,:));
subplot(4,1,3);plot(Time,X(3,:));
subplot(4,1,4);plot(Time,X(4,:));


%% EESA SERVO

A=[1 2 -1;1 0 1;-2 3 0];
B=[1 0;0 1;-1 -2];
M=[B A*B]
rank(M)
C=[1 0 0;0 1 0]
% y1=x1 y2=x2
% psi1=int(yr1-y1),  psi2=int(yr2-y2) psi=[psi1;psi2]
Ah=[A zeros(3,2);-C zeros(2)]
Bh=[B;zeros(2)]
M2=[Bh Ah*Bh Ah^2*Bh Ah^3*Bh]
rank(M2)
mu1=-5;mu2=-6;mu3=-7;mu4=-8;mu5=-9;
S1=null([Ah-mu1*eye(5) Bh])
S2=null([Ah-mu2*eye(5) Bh])
S3=null([Ah-mu3*eye(5) Bh])
S4=null([Ah-mu4*eye(5) Bh])
S5=null([Ah-mu5*eye(5) Bh])

v1q1=S1(:,2);v2q2=S2(:,1);v3q3=S3(:,2);v4q4=S4(:,2);v5q5=S5(:,2);
v1=v1q1(1:5);q1=v1q1(6:7);
v2=v2q2(1:5);q2=v2q2(6:7);
v3=v3q3(1:5);q3=v3q3(6:7);
v4=v4q4(1:5);q4=v4q4(6:7);
v5=v5q5(1:5);q5=v5q5(6:7);
K=-[q1 q2 q3 q4 q5]*inv([v1 v2 v3 v4 v5])
eig(Ah-Bh*K)


yr1 = @(t) 5*sign(sin(0.5*t));
yr2 = @(t) 10*sign(sin(0.25*t));
inp.K = K;
inp.A = A;
inp.B = B;
inp.C = C;
inp.y_ref1 = yr1;
inp.y_ref2 = yr2;
T=20;
dt=0.01;
Time = 0:dt:T;
X0=[0;0;0;0;0];

[X] = myRungeKutta1(@Fun2, Time, X0, inp);
figure;
subplot(3,1,1);plot(Time,X(1,:),'b',Time,yr1(Time),'r');
subplot(3,1,2);plot(Time,X(2,:),'b',Time,yr2(Time),'r');
subplot(3,1,3);plot(Time,X(3,:));




%% Functions

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
    yr1 = inp.y_ref1(t);
    yr2 = inp.y_ref2(t);
    y   = inp.C*X(1:end-2);
    u   = -inp.K*X;
    dX  = [inp.A*X(1:end-2)+inp.B*u ; yr1-y(1); yr2-y(2)];
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