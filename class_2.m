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
u0=1; % entrée, débit entrant en H1, converti en T/mn

%% Question 1: System matrices
A=[-1/T1 k2/T1
    k3/T2 -1/T2];
B=[(1-k1)/T1 k1/T2]';
C1=[1-k3 1-k2];        % Case 1
C2=[1 0];              % Case 2 
C3=[0 1];              % Case 3 
C4=[C1;C2];            % Case 4 
D = 0; % D4=[0; 0];

%% Question 2: 
% Stability analysis
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
step(Sys)

%%Question 3
%Controllability Test
