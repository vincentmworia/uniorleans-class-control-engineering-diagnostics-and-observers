clc
clf
close all 

%% Question 1: System matrices 

% system parameters
T1 = 5;
T2 = 1;
k1 = 0.5;
k2 = 0.3;
k3 = 0.1;
u0 = 1;

% System matrices
A = [-1/T1     k2/T1;
      k3/T2   -1/T2];

B = [(1-k1)/T1;
      k1/T2];

C1 = [1-k3  1-k2];   % Sensor 1: overall output
C2 = [1 0];          % Sensor 2: x1
C3 = [0 1];          % Sensor 3: x2
C4 = [C1; C2];       % Sensor 4 (not used here)

D = 0;

disp("System matrices:")
disp(A); disp(B); disp(C1); disp(C2); disp(C3);

n = size(A,1);   % number of states
disp("Number of states n = " + n);
 
%% Question 2: Stability analysis 

lambda = eig(A);

disp("Stability Test")
if all(real(lambda) < 0)
    disp("The system is stable")
else
    disp("The system is unstable")
end

% Transfer functions
G1 = tf(ss(A,B,C1,D));
G2 = tf(ss(A,B,C2,D));
G3 = tf(ss(A,B,C3,D));

disp("Transfer functions:")
disp(G1); disp(G2); disp(G3);

% Steady-state values
x_ss = -A\B*u0;
disp("Steady-state state values x_ss =")
disp(x_ss)

C = [C1; C2; C3];
y_ss = C*x_ss;
disp("Steady-state output values y_ss =")
disp(y_ss)

% Step response
Sys = ss(A,B,C,D);
figure
step(Sys)
title("Step response of the system")
 
%% Question 3: Controllability & Observability 

% ---- Controllability ----
disp("Controllability Test")
if rank(ctrb(A,B)) == n
    disp("The system is controllable")
else
    disp("The system is not controllable")
end
  
observability_test(A,C1,n,"Sensor 1 (C1)")
observability_test(A,C2,n,"Sensor 2 (C2)")
observability_test(A,C3,n,"Sensor 3 (C3)")

% Local function
function observability_test(A,C,n,name)

disp(newline +"Observability Test for " + name)

if rank(obsv(A,C)) == n
    disp("The system is observable")
else
    disp("The system is not observable")
end

end