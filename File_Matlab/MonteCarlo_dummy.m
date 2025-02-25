close all
set(0,'DefaultTextInterpreter','latex')
rng(1)

N      = 1e5; %Numero totale di particelle
eps    = 1e-2; %Parametro di riscalamento per il regime quasi-invariante
dt     = eps; %Passo di discretizzazione temporale
Tfin   = 1e1; %Tempo finale
nmax   = floor(Tfin/dt); %Numero di iterazioni temporali
lambda = 3.5;
sigma  = sqrt(6);
U      = -1+4*rand(N,1);
M10    = mean(U);

hbar = waitbar(0,'','Name','Iterazioni');
for n=1:nmax
    waitbar(n/nmax,hbar,sprintf('$n$ = %d / %d',n,nmax));
    U     = U(randperm(N));
    U1    = U(1:N/2); % per farle interagire
    U2    = U(N/2+1:end);
    Theta = binornd(1,dt/eps,N/2,1);
    eta   = sqrt(3)*(-1+2*rand(N/2,2)); 
    p     = 1-eps*lambda+sqrt(eps)*sigma*eta;
    q     = eps*lambda;
    U1new = (1-Theta).*U1+Theta.*(p(:,1).*U1+q*U2);
    U2new = (1-Theta).*U2+Theta.*(p(:,2).*U2+q*U1);
    U     = [U1new;
             U2new];
end
close(hbar)

uu   = min(U):0.01:max(U);
mu   = 2*lambda/sigma^2;
ginf = (mu*M10)^(1+mu)/gamma(1+mu).*exp(-mu*M10./uu)./(uu.^(2+mu));

figure(1)
histogram(U,'Normalization','pdf','LineStyle','-','FaceColor','#9ECB73','EdgeColor','#8CB665');
hold on
plot(uu,ginf,'LineWidth',1.9,'color',[0.9290, 0.6940, 0.1250]);
axis([-1.5 4.5 0 1.3]);
xlabel('$v$')
ylabel('$g^\infty(v)$')
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultLineLineWidth',1.2);
set(gca,'TickLabelInterpreter','latex')