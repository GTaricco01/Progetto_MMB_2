clearvars;close all;clc;
% contact matrix
load("covidItalia28Oct.mat");
M = [19.2 4.8 3.0 7.1 3.7 3.1 2.3 1.4 1.4;
     4.8 42.4 6.4 5.4 7.5 5.0 1.8 1.7 1.7;
     3.0 6.4 20.7 9.2 7.1 6.3 2.0 0.9 0.9;
     7.1 5.4 9.2 16.9 10.1 6.8 3.4 1.5 1.5;
     3.7 7.5 7.1 10.1 13.1 7.4 2.6 2.1 2.1;
     3.1 5.0 6.3 6.8 7.4 10.4 3.5 1.8 1.8;
     2.3 1.8 2.0 3.4 2.6 3.5 7.5 3.2 3.2;
     1.4 1.7 0.9 1.5 2.1 1.8 3.2 7.2 7.2;
     1.4 1.7 0.9 1.5 2.1 1.8 3.2 7.2 7.2];

infectiontime = find(confirmed);                   %indice dei giorni per cui si sono registrati contagi
timeVector    = timeVector(infectiontime);            %tolgo giorni per cui non avevo contagi in Italia
confirmed     = confirmed(infectiontime);              %tolgo i dati nulli
deaths        = deaths(infectiontime);
recovered     = recovered(infectiontime);
positive      = confirmed-deaths-recovered;  
removed       = recovered+deaths;

N_class = [441148; 5666380; 5962570; 6570438; 8061697; 9619515; 7964817; 6141543; 4543122];
N       = sum(N_class);       % Total population N = S + I + R
I0      = confirmed(1);       % initial number of infected (dai dati reali)
beta    = 8.45e-9;            % rate of infection
gamma   = 0.24;               % rate of recovery 
delta   = 0;                % rate of immunity loss
T       = length(timeVector); % same period of the real data [T] = [days]
% Inserire calcolo R0
R0      = beta*(N-I0)/gamma;
fprintf('Valore parametro R0 %.2f\n', R0)
X0   = [N,I0,0];
X0_c = zeros(9,3);
I0_c = zeros(9,1); I0_c(8) = I0;
% X0(:,1) = colonna dei suscettibili iniziali divisi per classe di età
% X0(:,2) = colonna degli infetti iniziali divisi per classe di età
% X0(:,3) = colonna dei rimossi iniziali divisi per classe di età

X0_c(:,1) = N_class;
X0_c(:,2) = I0_c;

% supponiamo di avere come infetti iniziali solo 2 individui della classe
% numero 8: per avere un'evoluzione significativa dell'epidemia nel tempo
% ho bisogno di implementare il contagio inter classe.


%%
p      = [beta, gamma];
% chiamo la funzione ode45 una volta per ogni classe
[t, X]     = ode45(@(t,x)sir_model(t,x,p,1,eye(9)),[0 T-1],X0(1,:));
% M(1:3,1:3) = M(1:3,1:3)./2;
[t_c, X_c] = ode45(@(t,x)sir_model(t,x,p,8,M),[0 T-1],X0_c(8,:));
S         = X(:,1);
I         = X(:,2);
R         = X(:,3);
S_c       = X_c(:,1);
I_c       = X_c(:,2);
R_c       = X_c(:,3);

figure(1)
subplot(1,2,1)
plot(t,I,'r','LineWidth',2);
xlabel('Days'); ylabel('Number of individuals'); title('Infected SIR model')
set(gca,'FontSize',10)
subplot(1,2,2)
plot(t_c,I_c,'r','LineWidth',2);
xlabel('Days'); ylabel('Number of individuals'); title('Infected SIR age model')
set(gca,'FontSize',10)

