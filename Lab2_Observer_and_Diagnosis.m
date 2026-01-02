% Master MARS
% Diagnosis and Observers
%
% Case study : Grinding-classification process

clc
%clf
close all
% system parameters
T1=5; % en mn
T2=1; % en mn
k1=0.5;
k2=0.3;
k3=0.1;
u0=1; % entrée, débit entrant en H1, converti en T/mn // Magnitude of input
n=2;
%% Question 1: System matrices

A=[-1/T1 k2/T1
    k3/T2 -1/T2];
B=[(1-k1)/T1 k1/T2]';
C1=[1-k3 1-k2];        % Case 1
C2=[1 0];              % Case 2 
C3=[0 1];              % Case 3 
C4=[C1;C2];            % Case 4 
D = 0; % D4=[0; 0];

%% Question 2: Stability analysis

lambda = eig(A);
if all(real(lambda)<0)
    display("The system is stable")
else
    display("The system is unstable")
end
%Transfer fuction

TF1=ss2tf(A,B,C1,D) % Case 1
TF2=ss2tf(A,B,C2,D) % Case 2 etc.

% Steady-state value :
x0 = -inv(A)*B*u0;

%Plotting the step response
C=[C1;C2; C3];
Sys=ss(A,B,C,D);
Sys1=ss(A,B,C1,D);
Sys2=ss(A,B,C2,D);
Sys3=ss(A,B,C3,D);

step(Sys)


%% Question 3: Stability analysis

%controlability test

Qc=ctrb(Sys)
rc=rank(Qc);
if rc==n
    display("The system is controllable")
else
    display("The system is not controllable")
end

%observability test
Qo=obsv(Sys)
ro=rank(Qo)

display("Observability Test for Output 1")
if ro==n
    display("The system is observable")
else
    display("The system is not observable")
end


%% Question 4: Luenberger Observer

%Case 1 (C1)
Cobs=C1;
% Desired Poles
p1=[-0.5 -1];
L = place(A',Cobs',p1); %Luenberger observer for sensor case1 (C1)
L=L'


%The observer as a system:
AL=A-L*Cobs; %observer state matrix
BL=[B L];  %input matrix (the observer has two inputs: u and y)
CL=Cobs;
sysobs=ss(AL, BL, CL, 0);

% System simulation to obtain the output y1 which is the input of the observer
[y,t,x] = step(Sys1);

%Observer Simulation
u=u0*ones(size(t1));

%Observer input
uobs=[u,y];

[yobs,tobs, xobs] = lsim(sysobs, uobs, t,[1;0.5]);
figure;
plot(tobs,xobs(:,1),'r',t,x(:,1),'b'); 
title('estimated and real value of the state: X1 (first component)')

figure;
plot(tobs,xobs(:,2),'r',t,x(:,2),'b'); 
title('estimated and real value of the state: X2 (second component)')

