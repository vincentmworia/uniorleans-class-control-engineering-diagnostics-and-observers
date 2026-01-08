% Master MARS
% Diagnostic et Observateurs
%
% Etude de cas : exemple du processus de broyage-classification
% One single script, grouped as Question 1 ... 8
%
% (Code kept exactly the same — only section headers added)

%% Question 1-3
% Master MARS
% Diagnostic et Observateurs
%
% Etude de cas : exemple du processus de broyage-classification
% Q. 1-3

%% Q.1 Modélisation
clc
%clf
close all
Broyeur0
% % Paramètres du système
% T1=5; % en mn
% T2=1; % en mn
% k1=0.5;
% k2=0.3;
% k3=0.1;
% u0=1/60; % entrée, débit entrant en H1, converti en T/mn
%
% A=[-1/T1 k2/T1
%     k3/T2 -1/T2];
% B=[(1-k1)/T1 k1/T2]';
% C1=[1-k3 1-k2];        % Config 1 des sorties
% C2=[1 0];              % Config 2 des sorties
% C3=[0 1];              % Config 3 des sorties
% C4=[C1;C2];            % Config 4 des sorties
% C=[C1;C2;C3];

%% Q.2. Étude de stabilité
disp('Valeurs propres de A')
lambda= eig(A)
if all(real(lambda)<=0),
    disp('Les valeurs propres de A sont à partie réelle négative, le sytème est donc stable')
else
    disp('A a au moins une valeur propre à partie réelle strictement positive, le sytème est donc instable')
end

%% Fonctions de transfert

%Systèmes définis dans Broyeur0 :
% sys1 = ss(A,B,C1,0);
% sys2= ss(A,B,C2,0);
% sys3= ss(A,B,C3,0);
% sys = ss(A,B,C,zeros(size(C,1),size(B,2))); %ss(A,B,C,0);
disp('Fonction de transfert config 1 (entre u et y):')
tf(sys1)
disp('Fonction de transfert config 2 (entre u et x1):')
tf(sys2)
disp('Fonction de transfert config 3 (entre u et x2):')
tf(sys3)
disp('Fonction de transfert toutes config (entre u et [y; x1; x2] :')
tf(sys)

%% Q.3 Commandabilité et observabilité
%Commandabilité :
n=2;
com = B;
Ocom=[];
for i=1:n-1,
    com = A*com;
    Ocom = [Ocom com]
end
%ou
Ocom_ =ctrb(sys)
rank(Ocom)
disp('Le rang de la matrice de commandabilité est 2 - le système est donc commandable')
%Observabilité config 1:
obs = C1;
Oobs=[];
for i=1:n-1,
    Oobs=[Oobs obs];
end
%ou
Oobs1=obsv(sys1)
rank(Oobs1)
disp('Le rang de la matrice d observabilité en config. 1 est 2 - le système est donc observable')
%Observabilité config 2:
obs = C2;
Oobs=[];
for i=1:n-1,
    Oobs=[Oobs obs];
end
%ou
Oobs2=obsv(sys2)
rank(Oobs2)
disp('Le rang de la matrice d observabilité en config. 2 est 2 - le système est donc observable')
%Observabilité config 3:
obs = C3;
Oobs=[];
for i=1:n-1,
    Oobs=[Oobs obs];
end
%ou
Oobs3=obsv(sys3)
rank(Oobs3)
disp('Le rang de la matrice d observabilité en config. 3 est 2 - le système est donc observable')

%Simulation
disp('Simulation la réponse indicielle')
t=0:0.2:30;
u=u0*ones(size(t));
y1=lsim(sys1,u,t);
disp('Tracé de y')
figure; lsim(sys1,u,t);pause
disp('Tracé de x1')
figure; lsim(sys2,u,t);pause
disp('Tracé de x2')
figure; lsim(sys3,u,t);pause


%% Question 4
% Master MARS
% Diagnostic et Observateurs
%
% Etude de cas : exemple du processus de broyage-classification
% QUESTION 4
%
% Reconstruction d'état à l'aide d'un observateur
%
disp('Reconstruction d etat par observateur de Luenberger')
%Données du système :
Broyeur0
%Défintion des paires de pôles :
p1=[-.5 -1]; % on désire des pôles à -1.5 et -2 pour l'observateur
p2=[-1.5 -2];
p3=[-2.5 -3];
%% pôles1 - config1
op=p1;
L=place(A',C1',op);
L=L';
%x0=[1;1]; %défini dans Broyeur0.m
%u=u0*ones(size(t));  %défini dans Broyeur0
%t=0:0.2:30; %défini dans Broyeur0
%Génération de y (une des entrées de l'observateur)
[y,t,x]=lsim(sys1,u',t,x0); % sys1 à définir aussi dans Broyeur0.m
%
% Définition de l'observateur en tant que système
AL=A-L*C1 %matrice d'état de l'observation
BL=[B L];
sysobs1p1=ss(AL,BL,C1,zeros(1,1));
Uobs=[u' y];
X0=[0;0]; % état initial de l'observateur
[ys,ts,xs] = lsim(sysobs1p1,Uobs,t,X0);
figure(1)
plot(t,ys,'r',t,y,'g'),title('Sortie globale config 1 pôles 1')
figure(2)
plot(t,x(:,1)-xs(:,1)),title('Erreur d''estimation 1 config 1 pôles 1')
figure(3)
plot(t,x(:,2)-xs(:,2)),title('Erreur d''estimation 2 config 1 pôles 1')

%% pôles2 - config1
op=p2;
L=place(A',C1',op);
L=L';
%x0=[1;1]; %défini dans Broyeur0.m
%u=u0*ones(size(t));  %défini dans Broyeur0
%t=0:0.2:30; %défini dans Broyeur0
%Génération de y (une des entrées de l'observateur)
[y,t,x]=lsim(sys1,u',t,x0); % sys1 à définir aussi dans Broyeur0.m
%
% Définition de l'observateur en tant que système
AL=A-L*C1 %matrice d'état de l'observation
BL=[B L];
sysobs1p1=ss(AL,BL,C1,zeros(1,1));
Uobs=[u' y];
X0=[0;0]; % état initial de l'observateur
[ys,ts,xs] = lsim(sysobs1p1,Uobs,t,X0);
figure;
plot(t,ys,'r',t,y,'g'),title('Sortie globale config 1 pôles 2')
figure;
plot(t,x(:,1)-xs(:,1)),title('Erreur d''estimation 1 config 1 pôles 2')
figure;
plot(t,x(:,2)-xs(:,2)),title('Erreur d''estimation 2 config 1 pôles 2')

%% pôles3 - config1
op=p3;
L=place(A',C1',op);
L=L';
%x0=[1;1]; %défini dans Broyeur0.m
%u=u0*ones(size(t));  %défini dans Broyeur0
%t=0:0.2:30; %défini dans Broyeur0
%Génération de y (une des entrées de l'observateur)
[y,t,x]=lsim(sys1,u',t,x0); % sys1 à définir aussi dans Broyeur0.m
%
% Définition de l'observateur en tant que système
AL=A-L*C1 %matrice d'état de l'observation
BL=[B L];
sysobs1p1=ss(AL,BL,C1,zeros(1,1));
Uobs=[u' y];
X0=[0;0]; % état initial de l'observateur
[ys,ts,xs] = lsim(sysobs1p1,Uobs,t,X0);
figure;
plot(t,ys,'r',t,y,'g'),title('Sortie globale config 1 pôles 3')
figure;
plot(t,x(:,1)-xs(:,1)),title('Erreur d''estimation 1 config 1 pôles 3')
figure;
plot(t,x(:,2)-xs(:,2)),title('Erreur d''estimation 2 config 1 pôles 3')


%% Question 5
% QUESTION 5
%Broyeur5

% Influence d'un défaut capteur
% Système de mesure complet
% Essayer avec les 2 configurations pour les pôles de l'observateur
% (-0.5, -1) et (-2.5, -3)
%
% Définition du système
%
Broyeur0
%u=u0*ones(size(t));
%t=0:0.2:30;
[y,t,x]=lsim(sys,u,t,x0);
%
% Simulation des défauts
%
y(30:33,1)=y(30:33,1)+0.5;  %5.8 - 6.4 s
y(60:63,2)=y(60:63,2)-0.25;  %11.8 - 12.4 s
%
% Reconstruction à l'aide d'un observateur : instrumentation C1
%
% Choix de pôles N° 1
op=[-0.5 -1] % on désire des pôles à -0.5 et -1  pour l'observateur

L=place(A',C1',op);
AL=A-L'*C1;
BL=[B L'];
%
% Définition de l'observateur en tant que système
%
y1=y(:,1);
sysobs1=ss(AL,BL,C1,zeros(1,2));
U=[u0*ones(size(t)) y1];
X0=[0;0]; % état initial
[ys,ts,xs] = lsim(sysobs1,U,t,X0);
figure(1)
plot(t,y),title('Sortie globale mesurée - poles1')
figure(2)
plot(t,xs(:,1),'r',t,xs(:,2),'b'),title('Etats estimés - poles1')
%
% Création d'un vecteur de sortie estimée
%
ye=xs*C';
%
figure(3)
plot(t,y(:,1)-ye(:,1)),title('Erreur de reconstruction 1 - poles1')
figure(4)
plot(t,y(:,2)-ye(:,2)),title('Erreur de reconstruction 2 - poles1')
figure(5)
plot(t,y(:,3)-ye(:,3)),title('Erreur de reconstruction 3 - poles1')

%
% Choix de pôles N° 2
op=[-2.5 -3] % on désire des pôles à -2.5 et -3 pour l'observateur

L=place(A',C1',op);
AL=A-L'*C1;
BL=[B L'];
%
% Définition de l'observateur en tant que système
%
y1=y(:,1);
sysobs1=ss(AL,BL,C1,zeros(1,2));
U=[u0*ones(size(t)) y1];
X0=[0;0]; % état initial
[ys,ts,xs] = lsim(sysobs1,U,t,X0);
figure(1)
plot(t,y),title('Sortie globale mesurée')
figure(2)
plot(t,xs(:,1),'r',t,xs(:,2),'b'),title('Etats estimés')
%
% Création d'un vecteur de sortie estimée
%
ye=xs*C';
%
figure(6)
plot(t,y(:,1)-ye(:,1)),title('Erreur de reconstruction 1')
figure(7)
plot(t,y(:,2)-ye(:,2)),title('Erreur de reconstruction 2')
figure(8)
plot(t,y(:,3)-ye(:,3)),title('Erreur de reconstruction 3')

% Fonctions de transfert des biais sur les erreurs de reconstruction

% op=[-2.5 -3]

L=place(A',C1',op);
AL=A-L'*C1;
BL=[B L'];
%
% Définition de l'observateur en tant que système
%
H=[1 0;0 1;0 0];
H1=H(1,:);
% Choix de pôles N° 1
op=[-0.5 -1] % on désire des pôles à -0.5 et -1 pour l'observateur
L=place(A',C1',op);
L=L'
Cdef=-C;
Adef=A-L*C1;
Bdef=L*H1;
Ddef=H;
[num1,den1]=ss2tf(Adef,Bdef,Cdef,Ddef,1)
[num2,den2]=ss2tf(Adef,Bdef,Cdef,Ddef,2)

% Choix de pôles N° 2
op=[-2.5 -3] % on désire des pôles à -2.5 et -3 pour l'observateur
L=place(A',C1',op);
L=L';
Cdef=-C;
Adef=A-L*C1;
Bdef=L*H1;
Ddef=H;
[num1,den1]=ss2tf(Adef,Bdef,Cdef,Ddef,1)
[num2,den2]=ss2tf(Adef,Bdef,Cdef,Ddef,2)


%% Question 6
% Question 6
%Broyeur6

% Influence du bruit de mesure
% Système de mesure C1
% Essayer avec les 2 configurations pour les pôles de l'observateur (-0.5, -1) et (-2.5, -3)
clearvars

Broyeur0
%
% définition du système
%
[y,t,x]=lsim(sys,u,t,x0);
%
% Simulation des défauts
%
y(30:33,1)=y(30:33,1)+0.5;
y(60:63,2)=y(60:63,2)-0.25;
%
% Simulation du bruit de mesure
%
b=randn(size(y)); % bruit normal
%b=(b-ones(length(y),1)*mean(b))./(ones(length(y),1)*std(b));
b=(b-mean(b))/std(b);
% à moyenne nulle et à variance unité
ym=y+0.15*b; % bruit de mesure de variance 0.15^2=0.0225
%
% Reconstruction à l'aide d'un observateur : instrumentation C1
%
op=[-10 -15]; % on désire des pôles à -0.5 (-2.5) et -1 (-3) pour l'observateur
L=place(A',C1',op);
L=L';
AL=A-L*C1;
BL=[B L];
%
% Définition de l'observateur en tant que système
%
y1=ym(:,1);
sysobs1=ss(AL,BL,C1,zeros(1,2));
U=[u0*ones(size(t)) y1];
X0=[0;0]; % état initial
[ys,ts,xs] = lsim(sysobs1,U,t,X0);
figure;
plot(t,ym),title('Sortie globale mesurée'); pause
figure;
plot(t,xs(:,1),'r',t,xs(:,2),'b'),title('Etats estimés'); pause
%
% Création d'un vecteur de sortie estimée
%
ye=xs*C';
figure;
plot(t,ym(:,1)-ye(:,1)),title('Erreur de reconstruction 1'); pause
figure;
plot(t,ym(:,2)-ye(:,2)),title('Erreur de reconstruction 2'); pause
figure;
plot(t,ym(:,3)-ye(:,3)),title('Erreur de reconstruction 3'); pause
%
% Recherche des fonctions de transfert entre y~ et b : (-C(sI-A+LC1)-1LH + I) où H=[1 0 0]
%
H=[1 0 0];
Ab=A-L*C1;
Bb=L*H;
Cb=-C;
Db=eye(3);
[num1,den1]=ss2tf(Ab,Bb,Cb,Db,1)
[num2,den2]=ss2tf(Ab,Bb,Cb,Db,2)
[num3,den3]=ss2tf(Ab,Bb,Cb,Db,3)
%
% Tracé du Bode
%
figure;
sysbr=tf(num1(1,:),den1)
bode(sysbr)


%% Question 7
%Broyeur 7
% Observateur multiple
% Système de mesure totale
%
% définition du système :
Broyeur0
[y,t,x]=lsim(sys,u,t,x0);
%
% Simulation des défauts
%
y(30:33,1)=y(30:33,1)+0.5;
y(60:63,2)=y(60:63,2)-0.25;
y(90:93,3)=y(90:93,2)+0.25;
%
H=eye(3);
H1=H(1,:);
H2=H(2,:);
H3=H(3,:);
H4=H(1:2,:);
%
%% Premier banc : instrumentation C1
disp('1er banc : instrumentation C1')
%
% Reconstruction à l'aide d'un observateur : instrumentation C1
%
op=[-1.5 -2]; % on désire des pôles à -1.5 et -2 pour Lobservateur
L=place(A',C1',op);
L=L';
AL=A-L*C1;
BL=[B L];
%
% Définition de Lobservateur en tant que système
%
y1=y(:,1);
sysobs1=ss(AL,BL,C1,zeros(1,2));
U=[u0*ones(size(t)) y1];
X0=[0;0]; % état initial
[ys,ts,xs] = lsim(sysobs1,U,t,X0);
figure(1)
subplot(511),plot(t,y),title('Sortie globale mesurée')
subplot(512),plot(t,xs(:,1),'r',t,xs(:,2),'b'),title('Etats estimés')
%
% Création d'un vecteur de sortie estimée
%
ye=xs*C';
subplot(513),plot(t,y(:,1)-ye(:,1)),title('Erreur de reconstruction 1')
subplot(514),plot(t,y(:,2)-ye(:,2)),title('Erreur de reconstruction 2')
subplot(515),plot(t,y(:,3)-ye(:,3)),title('Erreur de reconstruction 3')
%
%Transfert résidu / défaut :
%
Af=AL;
Bf=-L*H1;
Cf=C;
Df=H;
sysf=ss(Af,Bf,Cf,Df);
tf_f= tf(sysf)
pause
%
%% Deuxieme banc : instrumentation C2
disp('2eme banc : instrumentation C2')
%
% Reconstruction à Laide d'un observateur : instrumentation C2
%
op=[-1.5 -2]; % on désire des pôles à -1.5 et -2 pour Lobservateur
L=place(A',C2',op);
L=L';
AL=A-L*C2;
BL=[B L];
%
% Définition de Lobservateur en tant que système
%
y2=y(:,2);
sysobs1=ss(AL,BL,C2,zeros(1,2));
U=[u0*ones(size(t)) y2];
X0=[0;0]; % état initial
[ys,ts,xs] = lsim(sysobs1,U,t,X0);
figure(2)
subplot(511),plot(t,y),title('Sortie globale mesurée')
subplot(512),plot(t,xs(:,1),'r',t,xs(:,2),'b'),title('Etats estimés')
%
% Création d'un vecteur de sortie estimée
%
ye=xs*C';
subplot(513),plot(t,y(:,1)-ye(:,1)),title('Erreur de reconstruction 1')
subplot(514),plot(t,y(:,2)-ye(:,2)),title('Erreur de reconstruction 2')
subplot(515),plot(t,y(:,3)-ye(:,3)),title('Erreur de reconstruction 3')
%
%Transfert résidu / défaut :
%
Af=AL;
Bf=-L*H2;
Cf=C;
Df=H;
sysf=ss(Af,Bf,Cf,Df);
tf_f= tf(sysf)
pause
%
%% Troisieme banc : instrumentation C4
disp('3eme banc : instrumentation C4')
%
% Reconstruction à l'aide d'un observateur : instrumentation C4
%
op=[-1.5 -2]; % on désire des pôles à -1.5 et -2 pour Lobservateur
L=place(A',C4',op);
L=L';
AL=A-L*C4;
BL=[B L];
%
% Définition de Lobservateur en tant que système
%
y3=y(:,1:2);
sysobs1=ss(AL,BL,C4,zeros(2,3));
U=[u0*ones(size(t)) y3];
X0=[0;0]; % état initial
[ys,ts,xs] = lsim(sysobs1,U,t,X0);
figure(3)
subplot(511),plot(t,y),title('Sortie globale mesurée')
subplot(512),plot(t,xs(:,1),'r',t,xs(:,2),'b'),title('Etats estimés')
%
% Création d'un vecteur de sortie estimée
%
ye=xs*C';
subplot(513),plot(t,y(:,1)-ye(:,1)),title('Erreur de reconstruction 1')
subplot(514),plot(t,y(:,2)-ye(:,2)),title('Erreur de reconstruction 2')
subplot(515),plot(t,y(:,3)-ye(:,3)),title('Erreur de reconstruction 3')
%
%Transfert résidu / défaut :
%
Af=AL;
Bf=-L*H4;
Cf=C;
Df=H;
sysf=ss(Af,Bf,Cf,Df);
tf_f= tf(sysf)
pause
%
%% Quatrieme banc : instrumentation C3
disp('4eme banc : instrumentation C3')
%
% Reconstruction à Laide d'un observateur : instrumentation C3
%
op=[-1.5 -2]; % on désire des pôles à -1.5 et -2 pour Lobservateur
L=place(A',C3',op);
L=L';
AL=A-L*C3;
BL=[B L];
%
% Définition de Lobservateur en tant que système
%
y4=y(:,3);
sysobs1=ss(AL,BL,C3,zeros(1,2));
U=[u0*ones(size(t)) y4];
X0=[0;0]; % état initial
[ys,ts,xs] = lsim(sysobs1,U,t,X0);
figure(4)
subplot(511),plot(t,y),title('Sortie globale mesurée')
subplot(512),plot(t,xs(:,1),'r',t,xs(:,2),'b'),title('Etats estimés')
%
% Création d'un vecteur de sortie estimée
%
ye=xs*C';
subplot(513),plot(t,y(:,1)-ye(:,1)),title('Erreur de reconstruction 1')
subplot(514),plot(t,y(:,2)-ye(:,2)),title('Erreur de reconstruction 2')
subplot(515),plot(t,y(:,3)-ye(:,3)),title('Erreur de reconstruction 3')
%
%Transfert résidu / défaut :
%
Af=AL;
Bf=-L*H3;
Cf=C;
Df=H;
sysf=ss(Af,Bf,Cf,Df);
tf_f= tf(sysf)


%% Question 8
% Broyeur
% Influence d'une fuite
% Observateur à entrée inconnue
Broyeur0
%
% Définition du système
% Entrée du système avec fuite : [u f]
D=[0;-1]; % matrice d'action de la fuite généralisée
sys1f=ss(A,[B D],C,zeros(3,2));
%x0=[1;1];
%t=0:0.2:30;
%u=u0*ones(size(t));
%
% Création de la fuite
%
q=zeros(size(u));
tf1 = 8.8/0.2+1;
tf2=12.4/0.2+1;
q(tf1:tf2)=0.5;
f=(1/T2)*q+[diff(q) 0];
[y,t,x]=lsim(sys1f,[u;f],t,x0);
%
% Simulation des défauts
%
y(30:33,1)=y(30:33,1)+0.5;
y(60:63,2)=y(60:63,2)-0.25;
%
% Reconstruction à l'aide d'un observateur à entrée inconnue
%
N=[-0.5 0;0 -1]; % on désire des pôles à -0.5 et -1 pour l'observateur à entrée inconnue
CD_=pinv(C*D);
P=eye(2)-D*CD_*C;
G=P*B;
KC=P*A-N;
K=KC*pinv(C);
L=K+N*D*CD_;
%
% Définition de l'observateur en tant que système
%
sysobs1=ss(N,[G L],C,zeros(3,4));
U=[u0*ones(size(t)) y];
X0=[0;0]; % état initial
[ys,ts,z] = lsim(sysobs1,U,t,X0);
xs=(z'+D*pinv(C*D)*y')';
figure(1)
subplot(211),plot(t,y),title('Sortie globale mesurée')
subplot(212),plot(t,xs(:,1),'r',t,xs(:,2),'b'),title('Etats estimés')
%
% Création d'un vecteur de sortie estimée
%
ye=xs*C';
figure(2)
subplot(311),plot(t,y(:,1)-ye(:,1)),title('Erreur de reconstruction 1')
subplot(312),plot(t,y(:,2)-ye(:,2)),title('Erreur de reconstruction 2')
subplot(313),plot(t,y(:,3)-ye(:,3)),title('Erreur de reconstruction 3')
