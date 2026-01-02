clc; clf; close all;

% Parameters
T1 = 5; T2 = 1; k1 = 0.5; k2 = 0.3; k3 = 0.1;

% State-space (underflow = k; H1 overflow -> B1)
A = [-1/T1,   k2/T1;
      k3/T2, -1/T2];
B = [(1-k1)/T1;
      k1/T2];

% Outputs per case
C1 = [1-k3, 1-k2];  D1 = 0;   % Case 1: overall product
C2 = [1, 0];        D2 = 0;   % Case 2: x1
C3 = [0, 1];        D3 = 0;   % Case 3: x2

sys1 = ss(A,B,C1,D1);
sys2 = ss(A,B,C2,D2);
sys3 = ss(A,B,C3,D3);

%% (1) Stability: eigenvalues
disp('Eigenvalues of A:'); disp(eig(A));

%% (4) Transfer functions
disp('G1(s) (overall product):'); tf(sys1)
disp('G2(s) (x1):');             tf(sys2)
disp('G3(s) (x2):');             tf(sys3)

%% (5) Steady-state for unit step u0 = 1 (T/min)
u0 = 1;
x_ss = -inv(A)*B*u0;            % steady-state states
y1_ss = C1*x_ss + D1*u0;        % Case 1 output
y2_ss = C2*x_ss + D2*u0;        % Case 2 output
y3_ss = C3*x_ss + D3*u0;        % Case 3 output
disp('Steady-state x = [-A^{-1} B]*u0:'), disp(x_ss);
disp('Steady-state y (Case 1,2,3):'), disp([y1_ss, y2_ss, y3_ss]);

%% (6) Simulate unit step response
t = 0:0.1:50; u = ones(size(t));
[y1,~,x1] = lsim(sys1,u,t);
[y2,~,x2] = lsim(sys2,u,t);
[y3,~,x3] = lsim(sys3,u,t);

figure; plot(t,y1,'LineWidth',1.5); grid on; title('Case 1: Overall product y');
xlabel('time [min]'); ylabel('y [T/min]');

figure; plot(t,y2,'LineWidth',1.5); grid on; title('Case 2: x_1');
xlabel('time [min]'); ylabel('x_1 [T/min]');

figure; plot(t,y3,'LineWidth',1.5); grid on; title('Case 3: x_2');
xlabel('time [min]'); ylabel('x_2 [T/min]');

% (Optional) plot states for Case 1
figure; plot(t,x1,'LineWidth',1.5); grid on; title('States during Case 1 step');
legend('x_1','x_2'); xlabel('time [min]'); ylabel('flow [T/min]');