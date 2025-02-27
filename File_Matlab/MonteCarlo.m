function [f,U,n] = MonteCarlo(U0,beta,gamma,N)
set(0,'DefaultTextInterpreter','latex')
rng(1)

eps    = 1e-2; %Parametro di riscalamento per il regime quasi-invariante
dt     = eps; %Passo di discretizzazione temporale
Tfin   = 1e1; %Tempo finale
nmax   = floor(Tfin/dt); %Numero di iterazioni temporali

hbar = waitbar(0,'','Name','Iterazioni');
for n=1:nmax
    waitbar(n/nmax,hbar,sprintf('$n$ = %d / %d',n,nmax));
    % U     = U(randperm(N));
    U1    = U0(1:N/2); % per farle interagire le divido in due gruppi
    U2    = U0(N/2+1:end);
    Theta_b = binornd(1,beta,N/2,1);
    Theta_g = binornd(1,gamma,N/2,1);
    U1new = zeros(N/2,1);
    U2new = zeros(N/2,1);
    for p = 1:min(numel(U1),numel(U2))
        if (U1(p)<=-1) && (abs(U2(p))<=1) % interazione S-I
            UU1 = (1-beta)*U1(p) + beta*(U2(p)+2); % stati post interazione
            UU2 = U2(p);
            U1new(p) = (1-Theta_b(p)).*U1(p) + Theta_b(p).*(UU1); % aggiornamento temporale degli stati
            U2new(p) = (1-Theta_b(p)).*U2(p) + Theta_b(p).*(UU2);
        elseif (abs(U1(p))<=1) && (abs(U2(p))<=1) % interazione I-I
            UU1 = (1-gamma)*U1(p) + gamma*(U2(p)+2); % stati post interazione
            UU2 = U2(p);
            U1new(p) = (1-Theta_g(p)).*U1(p) + Theta_g(p).*(UU1);
            U2new(p) = (1-Theta_g(p)).*U2(p) + Theta_g(p).*(UU2);
        else
            U1new(p) = U1(p);
            U2new(p) = U2(p);
        end
    end
    U = [U1new;
         U2new];
end
close(hbar)

h = histogram(U,'LineStyle','-','FaceColor','#9ECB73','EdgeColor','#8CB665');
% ricavo la distribuzione prendendo i valori degli istogrammi
f = h.Values;
n = h.NumBins;
figure(1)
plot(h.BinLimits(1):h.BinWidth:h.BinLimits(2)-h.BinWidth,h.Values)
% legend(sprintf('Distribution of class i at time step %d',t))
xlabel('$u$')
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultLineLineWidth',1.2);
set(gca,'TickLabelInterpreter','latex')
end