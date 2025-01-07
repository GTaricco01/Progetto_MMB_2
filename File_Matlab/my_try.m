clearvars;close all;clc;
% Parameters
L  = 100;                          % Length of the spatial domain
Nx = 10000;                         % Number of spatial grid points
dx = L / (Nx - 1);                 % Spatial step size
x  = linspace(-L, L, Nx);          % Spatial grid
cour = 0.9;

c  = 1;                            % Transport speed
Nt = 5000;                          % Number of time steps
dt = cour*dx/abs(c);  
% T = 5;
% dt = T/5000;                     % Time step size
lam = dt/dx;

% Initial condition
u        = zeros(Nx,1);            % All individuals are susceptible, except some infects – otherwise infection will not start
u(x<=-1) = 1;

% Source term function S(u, x, t)
source = @(u, x, t) 0;

% Time loop
for n = 1:Nt
    u_new = u;
    
    % Upwind scheme for transport term
    for j = 2:Nx-1
        % ujp = u(j+1);
        % uj  = u(j);
        % ujm = u(j-1);
        % u_new(j) = u(j) - lam*(flusso(uj,ujp,abs(c),lam,4)-flusso(ujm,uj,abs(c),lam,4));
        u_new(j) = u(j) - c * dt / dx * (u(j) - u(j-1)) + dt * source(u(j), x(j), n * dt);
    end
    
    % new time step
    u = u_new;
    
    if mod(n, 50) == 0
        plot(x, u, 'b', 'LineWidth', 1.5);
        title(['Time = ', num2str(n*dt)]);
        xlabel('x');
        ylabel('u');
        drawnow;
        % pause
    end
end 

