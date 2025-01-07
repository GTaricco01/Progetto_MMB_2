clearvars;close all;clc

load covidItalia28Oct.mat

infectiontime = find(confirmed);                   %indice dei giorni per cui si sono registrati contagi

timeVector    = timeVector(infectiontime);            %tolgo giorni per cui non avevo contagi in Italia

confirmed     = confirmed(infectiontime);              %tolgo i dati nulli

deaths        = deaths(infectiontime);

recovered     = recovered(infectiontime);

positive      = confirmed-deaths-recovered;

removed       = recovered+deaths;



X_c = [441148; 5666380; 5962570; 6570438; 8061697; 9619515; 7964817; 6141543; 4543122];

N       = sum(X_c);       % Total population N = S + I + R

I0      = confirmed(1);       % initial number of infected (dai dati reali)

beta    = 8.45e-9;            % rate of infection

gamma   = 0.24;               % rate of recovery

delta   = 0;                % rate of immunity loss

T       = length(timeVector); % same period of the real data [T] = [days]

% Inserire calcolo R0
R0      = beta*(N-I0)/gamma;
fprintf('Valore parametro R0 %.2f\n', R0)
X0_c = zeros(9,3);
I0_c = zeros(9,1); I0_c(7) = I0;
% X0(:,1) = colonna dei suscettibili iniziali divisi per classe di età
% X0(:,2) = colonna degli infetti iniziali divisi per classe di età
% X0(:,3) = colonna dei rimossi iniziali divisi per classe di età
X0_c(:,1) = X_c;
X0_c(:,2) = I0_c;
M =  [19.2 4.8 3.0 7.1 3.7 3.1 2.3 1.4 1.4;
    4.8 42.4 6.4 5.4 7.5 5.0 1.8 1.7 1.7;
    3.0 6.4 20.7 9.2 7.1 6.3 2.0 0.9 0.9;
    7.1 5.4 9.2 16.9 10.1 6.8 3.4 1.5 1.5;
    3.7 7.5 7.1 10.1 13.1 7.4 2.6 2.1 2.1;
    3.1 5.0 6.3 6.8 7.4 10.4 3.5 1.8 1.8;
    2.3 1.8 2.0 3.4 2.6 3.5 7.5 3.2 3.2;
    1.4 1.7 0.9 1.5 2.1 1.8 3.2 7.2 7.2;
    1.4 1.7 0.9 1.5 2.1 1.8 3.2 7.2 7.2];

% elementi per la discretizzazione
L = 10;
n = 100;
dx = 2*L/n;
uu = linspace(-L,L,n)';
f = zeros(n,1); % distribuzione iniziale
f(uu<=-1) = 1/(L-1);
plot(uu,f)
TT = 5;
tt = linspace(0,TT);
dt = TT/100;
n_class = size(M,1);
X_c = X0_c;
for i = 1:n_class
    for j = 1:n_class
        for k = 3:n
            f(uu(k)) = f(uu(k))+dt*(M(i,j)*((beta+gamma)*f(uu(k-2))+(2-beta-gamma)*f(uu(k)))*X_c(j,2)+2*M(i,j)*(X_c(j,1)+X_c(j,3))*f(uu(k)));
        end
        plot(uu,f(uu))
    end
    X_c(i,1) = sum(f(uu<=-1));
    X_c(i,2) = sum(f(uu>=-1));
    X_c(i,3) = sum(f(uu>=1));
end
