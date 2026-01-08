% Diagnostics & Observers
% Vincent Mworia

%% Questions
% Q1: Correct A, B, C matrices ✅
% Q2: eigenvalues + stability + TFs + steady-state + step ✅
% Q3: controllability/observability tests ✅
% Q4: Luenberger observer gain via place + simulation ✅
% Q5: inject sensor bias faults + residual plots + compare 2 pole sets ✅
% Q6: add Gaussian noise + still detect bias via residuals ✅
% Q8: UIO design (P,N,G,K,L) + simulate + show residuals insensitive to leak ✅

clc; clf; close all;

%% Question 1: System matrices

% System parameters
T1 = 5;    T2 = 1;
k1 = 0.5;  k2 = 0.3;  k3 = 0.1;
u0 = 1;

% State-space matrices
A = [-1/T1    k2/T1;
      k3/T2  -1/T2];

B = [(1-k1)/T1;
      k1/T2];

% Output matrices (sensors)
C1 = [1-k3  1-k2];   % Sensor 1: overall output
C2 = [1 0];          % Sensor 2: x1
C3 = [0 1];          % Sensor 3: x2
C4 = [C1; C2];       % Sensor 4: [y1; x1]

D = 0;

disp("System matrices:")
display(A); display(B);
display(C1); display(C2); display(C3);

n = size(A,1);
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

y_ss = [C1; C2; C3]*x_ss;
disp("Steady-state output values y_ss =")
display(y_ss)

% Step response
figure;
step(ss(A,B,[C1;C2;C3],D));
title("Step response of the system")

%% Question 3: Controllability & Observability (short version)

disp("Controllability Test")
if rank(ctrb(A,B)) == n
    disp("The system is controllable")
else
    disp("The system is not controllable")
end

Cset   = {C1, C2, C3};
Clabel = {'Sensor 1 (C1)','Sensor 2 (C2)','Sensor 3 (C3)'};

for k = 1:3
    Cobs = Cset{k};
    disp(newline + "Observability Test for " + Clabel{k})
    if rank(obsv(A,Cobs)) == n
        disp("The system is observable")
    else
        disp("The system is not observable")
    end
end

%% Question 4: Luenberger Observer (short version)

p1 = [-0.5 -1];
xhat0 = [1; 0.5];

Cset   = {C1, C2, C3};
Clabel = {'C1','C2','C3'};

for k = 1:3
    Cobs = Cset{k};  label = Clabel{k};

    Sys1 = ss(A,B,Cobs,D);
    [y,t,x] = step(u0*Sys1);

    L = place(A',Cobs',p1)';                 % observer gain
    sysobs = ss(A-L*Cobs,[B L],eye(2),0);    % xhat output

    uobs = [u0*ones(size(t)) y];
    [~,~,xhat] = lsim(sysobs,uobs,t,xhat0);

    figure; plot(t,xhat(:,1),'r',t,x(:,1),'b'); grid on;
    legend('Estimated x1','Real x1','Location','best');
    title(['Observer (' label '): x1'])

    figure; plot(t,xhat(:,2),'r',t,x(:,2),'b'); grid on;
    legend('Estimated x2','Real x2','Location','best');
    title(['Observer (' label '): x2'])
end

%% Question 5: Sensor bias faults + Residuals (2 pole sets)

Call = [C1; C2; C3];
Cobs = C4;
t = linspace(0,30,2000)';   u = u0*ones(size(t));

bias1 = 0.5;   t1_on=5.8;  t1_off=6.4;
bias2 = -0.25; t2_on=11.8; t2_off=12.4;

[~,~,x] = lsim(ss(A,B,eye(2),0),u,t);
y_true = (Call*x')';

y_meas = (Cobs*x')';
m = zeros(size(y_meas));
m(t>=t1_on & t<=t1_off,1) = bias1;
m(t>=t2_on & t<=t2_off,2) = bias2;
y_faulty = y_meas + m;

xhat0 = [1;0.5];
uobs  = [u y_faulty];

for p = {[-1.5 -2.0], [-2.5 -3.0]}
    poles = p{1};
    L = place(A',Cobs',poles)';
    sysobs = ss(A-L*Cobs,[B L],eye(2),0);
    [~,~,xhat] = lsim(sysobs,uobs,t,xhat0);

    res = y_true - (Call*xhat')';

    figure; plot(t,res); grid on; hold on;
    xline(t1_on,'--'); xline(t1_off,'--');
    xline(t2_on,'--'); xline(t2_off,'--');
    title(sprintf('Residuals with faults - Poles{%.1f,%.1f}',poles(1),poles(2)));
    legend('r1 (overall y)','r2 (x1)','r3 (x2)','Location','best');
end

%% Question 6: Measurement noise (Gaussian) + Residuals (2 pole sets)

Call = [C1; C2; C3];
Cobs = Call;
t = linspace(0,30,2000)';   u = u0*ones(size(t));

[~,~,x] = lsim(ss(A,B,eye(2),0),u,t);
y_true = (Call*x')';

m = zeros(size(y_true));
m(t>=t1_on & t<=t1_off,1) = bias1;
m(t>=t2_on & t<=t2_off,2) = bias2;

rng(1);
y_meas = y_true + m + sqrt(0.0225)*randn(size(y_true));

xhat0 = [1;0.5];
for p = {[-1.5 -2.0], [-2.5 -3.0]}
    poles = p{1};
    L = place(A',Cobs',poles)';
    sysobs = ss(A-L*Cobs,[B L],eye(2),0);
    [~,~,xhat] = lsim(sysobs,[u y_meas],t,xhat0);

    res = y_true - (Call*xhat')';

    figure; plot(t,res); grid on; hold on;
    xline(t1_on,'--'); xline(t1_off,'--');
    xline(t2_on,'--'); xline(t2_off,'--');
    title(sprintf('Residuals with faults + noise - Poles{%.1f,%.1f}',poles(1),poles(2)));
    legend('r1 (overall y)','r2 (x1)','r3 (x2)','Location','best');
end

%% Question 8: Unknown Input Observer (UIO) for leak (FINAL, exam-safe)

t = linspace(0,30,2000)';   u = u0*ones(size(t));

% Leak (unknown input): q(t)=0.5 between 8.8 and 12.4  -> f(t)=q/T2 (dq/dt=0)
tL_on=8.8; tL_off=12.4; q0=0.5;
f = (q0/T2) * (t>=tL_on & t<=tL_off);

Cui = [C1; C2; C3];     % measured outputs used for UIO
Dui = [0; -1];          % leak distribution (acts on x2 dynamics as loss)

% --- UIO existence condition check: rank(CD)=rank(D) ---
if rank(Cui*Dui) ~= rank(Dui)
    error("UIO condition fails: rank(C*D) must equal rank(D). Choose another sensor set.");
end

% --- Plant with unknown input (u and f) ---
Sysxf = ss(A,[B Dui],eye(2),0);
[~,~,x] = lsim(Sysxf,[u f],t);
y = (Cui*x')';          % true outputs (contain leak through plant)

% Add sensor bias faults (same as Q5) to measurements (faults only, leak is NOT a fault here)
mf = zeros(size(y));
mf(t>=t1_on & t<=t1_off,1) = bias1;   % sensor 1 bias
mf(t>=t2_on & t<=t2_off,2) = bias2;   % sensor 2 bias
y_meas = y + mf;

% --- UIO design (prof-style formulas, robust + standard) ---
N    = diag([-1.5 -2.0]);             % stable observer dynamics
CD_  = pinv(Cui*Dui);
P    = eye(2) - Dui*CD_*Cui;
G    = P*B;
K    = (P*A - N)*pinv(Cui);
Lui  = K + N*Dui*CD_;

% --- UIO simulation: zdot = N z + G u + Lui y,  xhat = z + D(CD)^+ y ---
sysZ = ss(N,[G Lui],eye(2),0);        % inputs: [u ; y1 y2 y3]
[~,~,z] = lsim(sysZ,[u y_meas],t,[0;0]);

xhat = z + y_meas*(Dui*CD_)';         % estimated states
yhat = (Cui*xhat')';                  % estimated outputs
r    = y_meas - yhat;                 % residuals

figure; plot(t,r); grid on; hold on;
xline(t1_on,'--'); xline(t1_off,'--');
xline(t2_on,'--'); xline(t2_off,'--');
xline(tL_on,':');  xline(tL_off,':');
title('UIO residuals: sensor faults visible, leak rejected');
legend('r1 (overall y)','r2 (x1)','r3 (x2)','Location','best');
