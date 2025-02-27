clearvars; close all;clc
% particelle per MonteCarlo
N_class = [441148; 5666380; 5962570; 6570438; 8061697; 9619515; 7964817; 6141543; 4543122];
% N = sum(N_class);
N = 1e5;
M =  [19.2 4.8 3.0 7.1 3.7 3.1 2.3 1.4 1.4;
       4.8 42.4 6.4 5.4 7.5 5.0 1.8 1.7 1.7;
       3.0 6.4 20.7 9.2 7.1 6.3 2.0 0.9 0.9;
       7.1 5.4 9.2 16.9 10.1 6.8 3.4 1.5 1.5;
       3.7 7.5 7.1 10.1 13.1 7.4 2.6 2.1 2.1;
       3.1 5.0 6.3 6.8 7.4 10.4 3.5 1.8 1.8;
       2.3 1.8 2.0 3.4 2.6 3.5 7.5 3.2 3.2;
       1.4 1.7 0.9 1.5 2.1 1.8 3.2 7.2 7.2;
       1.4 1.7 0.9 1.5 2.1 1.8 3.2 7.2 7.2];

L = 10;
Nu   = 1000;                            % Numero di punti spaziali
Tmax = 50;                            % Tempo massimo di simulazione
CFL  = 0.9;                            % Numero di Courant-Friedrichs-Lewyc
w    = @(u) 1;
% condizione iniziale
U0      = -10+9.1*rand(N,1);
% parametri
beta   = .5;  
gamma  = 0.24;
for t = 1:100
    % for N_c = N_class
        % primo passo
        [f_new_tilda,U,n,edges] = MonteCarlo(U0,beta,gamma,N);
        % secondo passo
        f_new = PassoUpwind(L,n,Tmax,CFL,w,f_new_tilda);
        % aggiornamento
        U0 = U;
        S = sum(f_new(1:find(edges==-1)));
        I = sum(f_new(find(edges==-1):find(edges==1)));
        R = sum(f_new(find(edges==1):end));
        sprintf('Suscettibili al tempo %d: %f. Infetti al tempo %d: %f. Rimossi al tempo %d: %f.',t,S,t,I,t,R)
    % end
end