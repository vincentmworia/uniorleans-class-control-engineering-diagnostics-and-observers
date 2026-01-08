T1=5; T2=1;
k1=0.5; k2=0.3; k3=0.1;
u0=1;                 % (use u0=1/60 if prof uses minutes)

A=[-1/T1 k2/T1
    k3/T2 -1/T2];

B=[(1-k1)/T1
    k1/T2];

C1=[1-k3 1-k2];
C2=[1 0];
C3=[0 1];
C4=[C1;C2];
C=[C1;C2;C3];

sys1=ss(A,B,C1,0);
sys2=ss(A,B,C2,0);
sys3=ss(A,B,C3,0);
sys =ss(A,B,C,zeros(3,1));

t=0:0.2:30;
u=u0*ones(size(t));
x0=[1;1];
