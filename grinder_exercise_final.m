% Diagnostics & Observers
% Vincent Mworia

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

y_ss = [C1; C2; C3]*x_ss;
disp("Steady-state output values y_ss =")
display(y_ss)

% Step response 
figure;
step(ss(A,B,C,D));
title("Step response of the system")
 
%% Question 3: Controllability & Observability (short version)

% Controllability
disp("Controllability Test")
if rank(ctrb(A,B)) == n
    disp("The system is controllable")
else
    disp("The system is not controllable")
end

% Observability
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

p1 = [-0.5 -1];            % desired observer poles
xhat0 = [1; 0.5];          % observer initial condition

Cset   = {C1, C2, C3};
Clabel = {'C1','C2','C3'};

for k = 1:3
    Cobs = Cset{k};  label = Clabel{k};

    % Real system for this sensor
    Sys1 = ss(A,B,Cobs,D);
    [y,t,x] = step(u0*Sys1);

    % Observer gain and observer system
    L = place(A',Cobs',p1)';                 % observer gain
    sysobs = ss(A-L*Cobs,[B L],eye(2),0);    % xhat output

    % Observer simulation (inputs: u and y)
    uobs = [u0*ones(size(t)) y];
    [~,~,xhat] = lsim(sysobs,uobs,t,xhat0);

    % Plots (x1 and x2)
    figure; plot(t,xhat(:,1),'r',t,x(:,1),'b'); grid on;
    legend('Estimated x1','Real x1','Location','best');
    title(['Observer (' label '): x1'])

    figure; plot(t,xhat(:,2),'r',t,x(:,2),'b'); grid on;
    legend('Estimated x2','Real x2','Location','best');
    title(['Observer (' label '): x2'])
end


%% Question 5: Sensor bias faults + Residuals (2 pole sets)

Call = [C1; C2; C3];        % outputs we want residuals for (3 outputs)
Cobs = C4;                  % measured outputs used by observer (y1 and y2)
t = linspace(0,30,2000)';   u = u0*ones(size(t));

% Fault definition (bias)
bias1 = 0.5;   t1_on=5.8;  t1_off=6.4;     % sensor 1 (y1)
bias2 = -0.25; t2_on=11.8; t2_off=12.4;    % sensor 2 (y2)

% Plant simulation -> get x(t)
[~,~,x] = lsim(ss(A,B,eye(2),0),u,t);

% True outputs (3 outputs)
y_true = (Call*x')';

% Measured outputs (2 outputs) + add faults
y_meas = (Cobs*x')';
m = zeros(size(y_meas));
m(t>=t1_on & t<=t1_off,1) = bias1;
m(t>=t2_on & t<=t2_off,2) = bias2;
y_faulty = y_meas + m;

% Observer initial condition + inputs
xhat0 = [1;0.5];
uobs  = [u y_faulty];

% Run observer for both pole sets
for p = {[-1.5 -2.0], [-2.5 -3.0]}
    poles = p{1};
    L = place(A',Cobs',poles)';               % observer gain
    sysobs = ss(A-L*Cobs,[B L],eye(2),0);     % xhat output
    [~,~,xhat] = lsim(sysobs,uobs,t,xhat0);

    res = y_true - (Call*xhat')';             % residuals

    figure; plot(t,res); grid on; hold on;
    xline(t1_on,'--'); xline(t1_off,'--');
    xline(t2_on,'--'); xline(t2_off,'--');
    title(sprintf('Residuals with faults - Poles{%.1f,%.1f}',poles(1),poles(2)));
    legend('r1 (overall y)','r2 (x1)','r3 (x2)','Location','best');
end



%% Question 6: Measurement noise (Gaussian) + Residuals (2 pole sets
Call = [C1; C2; C3];        % outputs to plot residuals for
Cobs = Call;                % sensors 1,2,3 are measured
t = linspace(0,30,2000)';   u = u0*ones(size(t));

% Simulate plant to get x(t) and true outputs
[~,~,x] = lsim(ss(A,B,eye(2),0),u,t);
y_true = (Call*x')';

% Add same bias faults as Q5 (on outputs 1 and 2)
m = zeros(size(y_true));
m(t>=t1_on & t<=t1_off,1) = bias1;
m(t>=t2_on & t<=t2_off,2) = bias2;

% Add Gaussian noise b(t) ~ N(0,0.0225)
rng(1);
y_meas = y_true + m + sqrt(0.0225)*randn(size(y_true));

% Run observer for both pole sets
xhat0 = [1;0.5];
for p = {[-1.5 -2.0], [-2.5 -3.0]}
    poles = p{1};
    L = place(A',Cobs',poles)';                 % observer gain
    sysobs = ss(A-L*Cobs,[B L],eye(2),0);       % xhat output
    [~,~,xhat] = lsim(sysobs,[u y_meas],t,xhat0);

    res = y_true - (Call*xhat')';               % residuals

    figure; plot(t,res); grid on; hold on;
    xline(t1_on,'--'); xline(t1_off,'--');
    xline(t2_on,'--'); xline(t2_off,'--');
    title(sprintf('Residuals with faults + noise - Poles{%.1f,%.1f}',poles(1),poles(2)));
    legend('r1 (overall y)','r2 (x1)','r3 (x2)','Location','best');
end

%% Question 8: UIO for leak (unknown input) + show residuals insensitive to leak

t = linspace(0,30,2000)';   u = u0*ones(size(t));

% Leak: q(t)=0.5 between 8.8 and 12.4  -> f(t)=q/T2 (dq/dt=0)
tL_on=8.8; tL_off=12.4; q0=0.5;
f = q0/T2 * (t>=tL_on & t<=tL_off);

C = [C1; C2; C3];           % 3 measured outputs
F = [0; -1];                % leak enters x2 dynamics (loss)

% --- Plant with unknown input (u and f) ---
Sysxf = ss(A,[B F],eye(2),0);
[~,~,x] = lsim(Sysxf,[u f],t);
y = (C*x')';

% Add sensor bias faults (same as Q5) to measurements
m = zeros(size(y));
m(t>=t1_on & t<=t1_off,1) = bias1;
m(t>=t2_on & t<=t2_off,2) = bias2;
y_meas = y + m;

% --- UIO design (course formulas) ---
N = diag([-1.5 -2.0]);          % chosen stable N
E = -F*pinv(C*F);               % E = -F(CF)^+
P = eye(2) + E*C;               % P = I + EC
G = P*B;                        % G = PB
K = (P*A - N)*pinv(C);          % K = (PA-N)C^+
L = K + N*F*pinv(C*F);          % L = K + N F (CF)^+

% --- UIO simulation: zdot = Nz + Gu + Ly,  xhat = z - Ey ---
sysZ = ss(N,[G L],eye(2),0);    % inputs: [u ; y1 y2 y3]
[~,~,z] = lsim(sysZ,[u y_meas],t,xhat0);
xhat = z - y_meas*E';

% Residuals: r = y_meas - yhat
yhat = (C*xhat')';
r = y_meas - yhat;

figure; plot(t,r); grid on; hold on;
xline(t1_on,'--'); xline(t1_off,'--');
xline(t2_on,'--'); xline(t2_off,'--');
xline(tL_on,':');  xline(tL_off,':');
title('UIO residuals: sensor faults visible, leak rejected');
legend('r1 (overall y)','r2 (x1)','r3 (x2)','Location','best');
