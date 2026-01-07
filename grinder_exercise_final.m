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
display(A); display(B); display(C1); display(C2); display(C3);

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
display(G1); display(G2); display(G3);

% Steady-state values
x_ss = -A\B*u0;
disp("Steady-state state values x_ss =")
display(x_ss)

C = [C1; C2; C3];
y_ss = C*x_ss;
disp("Steady-state output values y_ss =")
display(y_ss)

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

%% Question 4: Luenberger Observer 
 
p1 = [-0.5 -1];      % desired observer poles
xhat0 = [1; 0.5];    % observer initial condition

run_observer_case(A,B,C1,D,u0,p1,xhat0,"C1");
run_observer_case(A,B,C2,D,u0,p1,xhat0,"C2");
run_observer_case(A,B,C3,D,u0,p1,xhat0,"C3"); 

function run_observer_case(A,B,Cobs,D,u0,p1,xhat0,label)

% Observer gain
L = place(A', Cobs', p1);
L = L';
disp("Observer gain L (" + label + ") = ")
disp(L)

% Observer system
AL = A - L*Cobs;
BL = [B L];
CL = Cobs;

sysobs = ss(AL, BL, CL, 0);

% Real system
Sys1 = ss(A, B, Cobs, D);

% Generate output and true states
[y, t, x] = step(u0*Sys1);

% Observer input
u = u0*ones(size(t));
uobs = [u y];

% Observer simulation
[~, tobs, xobs] = lsim(sysobs, uobs, t, xhat0);

% Plots
figure;
plot(tobs, xobs(:,1),'r', t, x(:,1),'b');
grid on;
legend('Estimated x_1','Real x_1','Location','best');
title("Observer (" + label + "): x_1")

figure;
plot(tobs, xobs(:,2),'r', t, x(:,2),'b');
grid on;
legend('Estimated x_2','Real x_2','Location','best');
title("Observer (" + label + "): x_2")

end


%% Question 5: Influence of a sensor fault (bias) + Residuals 
Call = [C1; C2; C3];     % outputs we want residuals for (3 outputs)
Cobs = C4;               % measured outputs used by observer (2 outputs: y1 and y2)

% Fault definition (bias)
t = linspace(0,30,2000)';          
u = u0*ones(size(t));

bias1 = 0.5;   t1_on=5.8;  t1_off=6.4;    % sensor 1
bias2 = -0.25; t2_on=11.8; t2_off=12.4;   % sensor 2

% Plant simulation -> get x(t)
Sysx = ss(A,B,eye(2),0);
[~,~,x] = lsim(Sysx,u,t);

% True outputs (3 outputs)
y_true = (Call*x')';

% Measured outputs (2 outputs) + add faults
y_meas = (Cobs*x')';
m = zeros(size(y_meas));
m(t>=t1_on & t<=t1_off,1) = bias1;
m(t>=t2_on & t<=t2_off,2) = bias2;
y_faulty = y_meas + m;

xhat0 = [1;0.5];

%% --- Case 1 poles {-1.5,-2.0} ---
p = [-1.5 -2.0];
L = place(A',Cobs',p)';                  % observer gain
sysobs = ss(A-L*Cobs,[B L],eye(2),0);    % observer system
uobs = [u y_faulty];
[~,~,xhat] = lsim(sysobs,uobs,t,xhat0);

y_hat = (Call*xhat')';
res = y_true - y_hat;

figure;
plot(t,res(:,1)); hold on; grid on;
plot(t,res(:,2));
plot(t,res(:,3));
xline(t1_on,'--'); xline(t1_off,'--');
xline(t2_on,'--'); xline(t2_off,'--');
title('Residuals with faults - Poles{-1.5,-2.0}');
legend('r1 (overall y)','r2 (x1)','r3 (x2)','Location','best');

%% --- Case 2 poles {-2.5,-3.0} ---
p = [-2.5 -3.0];
L = place(A',Cobs',p)';
sysobs = ss(A-L*Cobs,[B L],eye(2),0);
[~,~,xhat] = lsim(sysobs,uobs,t,xhat0);

y_hat = (Call*xhat')';
res = y_true - y_hat;

figure;
plot(t,res(:,1)); hold on; grid on;
plot(t,res(:,2));
plot(t,res(:,3));
xline(t1_on,'--'); xline(t1_off,'--');
xline(t2_on,'--'); xline(t2_off,'--');
title('Residuals with faults - Poles{-2.5,-3.0}');
legend('r1 (overall y)','r2 (x1)','r3 (x2)','Location','best');


%% Question 6: Influence of measurement noise (Gaussian) + Residuals

% We keep the same faults as Q5 (bias on sensor 1 and 2),
% and now add Gaussian measurement noise b(t) to sensors 1,2,3.

Call = [C1; C2; C3];     % 3 outputs for residuals: [y1; y2; y3]
Cobs = Call;             % now we assume sensors 1,2,3 are measured (as in statement)

% Noise parameters
var_b = 0.0225;
sigma_b = sqrt(var_b);

% Use SAME simulation time and input
t = linspace(0,30,2000)';          
u = u0*ones(size(t));

% --- Plant simulation -> get x(t)
Sysx = ss(A,B,eye(2),0);
[~,~,x] = lsim(Sysx,u,t);

% True outputs (3 outputs)
y_true = (Call*x')';

% --- Build sensor faults m(t) (only on sensor 1 and 2)
m = zeros(size(y_true));                 % Nx3
m(t>=t1_on & t<=t1_off,1) = bias1;       % bias on output 1
m(t>=t2_on & t<=t2_off,2) = bias2;       % bias on output 2

% --- Add measurement noise b(t) to all 3 sensors
rng(1);                                  % for repeatable results (optional)
b = sigma_b*randn(size(y_true));         % Nx3 Gaussian noise

% Measured outputs with faults + noise
y_meas = y_true + m + b;

% Observer initial condition
xhat0 = [1;0.5];

%% --- Case 1 poles {-1.5,-2.0} ---
p = [-1.5 -2.0];
L = place(A',Cobs',p)';                    % observer gain

% Observer system: xhatdot = (A-LC)xhat + Bu + L*y_meas
sysobs = ss(A-L*Cobs,[B L],eye(2),0);

uobs = [u y_meas];
[~,~,xhat] = lsim(sysobs,uobs,t,xhat0);

y_hat = (Call*xhat')';
res = y_true - y_hat;

figure;
plot(t,res(:,1)); hold on; grid on;
plot(t,res(:,2));
plot(t,res(:,3));
xline(t1_on,'--'); xline(t1_off,'--');
xline(t2_on,'--'); xline(t2_off,'--');
title('Residuals with faults + noise - Poles{-1.5,-2.0}');
legend('r1 (overall y)','r2 (x1)','r3 (x2)','Location','best');

%% --- Case 2 poles {-2.5,-3.0} ---
p = [-2.5 -3.0];
L = place(A',Cobs',p)';

sysobs = ss(A-L*Cobs,[B L],eye(2),0);

[~,~,xhat] = lsim(sysobs,uobs,t,xhat0);

y_hat = (Call*xhat')';
res = y_true - y_hat;

figure;
plot(t,res(:,1)); hold on; grid on;
plot(t,res(:,2));
plot(t,res(:,3));
xline(t1_on,'--'); xline(t1_off,'--');
xline(t2_on,'--'); xline(t2_off,'--');
title('Residuals with faults + noise - Poles{-2.5,-3.0}');
legend('r1 (overall y)','r2 (x1)','r3 (x2)','Location','best');
