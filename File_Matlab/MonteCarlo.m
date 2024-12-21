close all
set(0,'DefaultTextInterpreter','latex')

N      = 1e6; %Numero totale di particelle
eps    = 1e-2; %Parametro di riscalamento per il regime quasi-invariante
dt     = eps; %Passo di discretizzazione temporale
Tfin   = 1e1; %Tempo finale
nmax   = floor(Tfin/dt); %Numero di iterazioni temporali
lambda = 3.5;
sigma  = sqrt(6);
V      = -1+4*rand(N,1);
M10    = mean(V);

hbar=waitbar(0,'','Name','Iterazioni');
for n=1:nmax
    waitbar(n/nmax,hbar,sprintf('$n$ = %d / %d',n,nmax));
    V     = V(randperm(N));
    V1    = V(1:N/2);
    V2    = V(N/2+1:end);
    Theta = binornd(1,dt/eps,N/2,1);
    eta   = sqrt(3)*(-1+2*rand(N/2,2)); 
    p     = 1-eps*lambda+sqrt(eps)*sigma*eta;
    q     = eps*lambda;
    V1new = (1-Theta).*V1+Theta.*(p(:,1).*V1+q*V2);
    V2new = (1-Theta).*V2+Theta.*(p(:,2).*V2+q*V1);
    V     = [V1new;
             V2new];
end
close(hbar)

vv   = min(V):0.01:max(V);
mu   = 2*lambda/sigma^2;
ginf = (mu*M10)^(1+mu)/gamma(1+mu).*exp(-mu*M10./vv)./(vv.^(2+mu));

figure(1)
histogram(V,'Normalization','pdf','LineStyle','-','FaceColor','#9ECB73','EdgeColor','#8CB665');
hold on
plot(vv,ginf,'LineWidth',1.9,'color',[0.9290, 0.6940, 0.1250]);
axis([-1.5 4.5 0 1.2717]);
xlabel('$v$')
ylabel('$g^\infty(v)$')
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultLineLineWidth',1.2);
set(gca,'TickLabelInterpreter','latex')