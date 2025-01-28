% questo è il caso nonlineare quindi:
% u_t + f(u)_x = 0
% u(x,0)      = u0(x)

clearvars, close all;

% parametri
T    = 2;
x_sx = -10;
x_dx = 10;

f  = @(x) 0.5.* x.^2;  % Burgers (flusso)
df = @(x) x;

%  f  = @(x) x.*(1-x);  % traffic
%  df = @(x) -2*x;

% griglie
dx = 0.01;
dt = 0.002;
X  = x_sx:dx:x_dx-dx;
tt = 0:dt:T;
a  = @(u,v,f) .5;

% u0=@(x) sin(2*pi*x);
u0=@(x) 1*((x>-10)&(x<0.1));
% u0=@(x) (-1+5*x).*(x>0.2&x<=0.4)+(x>0.4&x<=0.6)+(4-5*x).*(x>0.6&x<=0.8);

% implementazione con condizioni perdiodiche

U1=zeros(length(X),length(tt));
U1(:,1)=u0(X');

U2=zeros(length(X),length(tt));
U2(:,1)=u0(X');

vib1=@(X,Y,dt) X-sqrt(2*dt*mu(X))-Y;
for t=1:length(tt)-1 % ciclo sui tempi
    % fa il ciclo sulle celle on-line, cioè su una sola riga
    % Upwind
    U1(:,t+1) = U1(:,t)-0.5*dt/dx*( f([U1(2:end,t); U1(1,t)])-f([U1(end,t); U1(1:end-1,t)]))...
        +dt/dx*0.5*(abs(a(U1(:,t),[U1(2:end,t); U1(1,t)],f)).*([U1(2:end,t); U1(1,t)]-U1(:,t))...
        -abs(a([U1(end,t); U1(1:end-1,t)],U1(:,t),f)).*(U1(:,t)-[U1(end,t); U1(1:end-1,t)]))+dt/dx*.5;
        % le ultime due righe servono da stabilizzazione? 'a' è la velocità
        % di trasporto da definire come w(u) = w_tilda*H(1-u^2)

    % Lax-Wendroff
    U2(:,t+1)=U2(:,t)-0.5*dt/dx*( f([U2(2:end,t); U2(1,t)])-f([U2(end,t); U2(1:end-1,t)]))...
        +0.5*(dt/dx)^2*( (a(U2(:,t),[U2(2:end,t); U2(1,t)],f)).^2.*([U2(2:end,t); U2(1,t)]-U2(:,t)) ...
        -(a([U2(end,t); U2(1:end-1,t)],U2(:,t),f)).^2.*(U2(:,t)-[U2(end,t); U2(1:end-1,t)]))+dt/dx*.5;


    CFL=max(abs(a(U2(:,t),[U2(2:end,t); U2(1,t)],f)))*dt/dx;
    plot(X,U1(:,t+1),'o-b',X,U2(:,t+1),'+-g');
    legend('Upwind','Lax-Wendroff')
    % axis([x_sx x_dx min(U2(:,1))-0.3 max(U2(:,1))+0.3])
    title(['CFL =' num2str(CFL)]);
    pause(0.01)
end





