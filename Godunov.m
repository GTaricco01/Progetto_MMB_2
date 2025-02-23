clear; clc; close all;

% Parametri
L = 10;                               % Lunghezza del dominio
Nx = 1000;                            % Numero di punti spaziali
Tmax = 10;                            % Tempo massimo di simulazione
CFL = 0.9;                            % Numero di Courant-Friedrichs-Lewy
w = -0.5;                             % Velocità di drift (negativa)

% Discretizzazione spaziale
dx = 2*L / Nx;
x = linspace(-L, L, Nx);


% Condizione iniziale: onda quadrata
f0 = @(u) 2 * (abs(u)<1)+ 10*(u<=-1); 
u = zeros(1, Nx);

for i = 1:Nx-1
    u(i) = (f0(x(i)) + f0(x(i+1))) / 2; % Media per inizializzazione
end

% Calcolo del passo temporale (Condizione CFL corretta)
dt = CFL * (dx /abs(w));
Nt = ceil(Tmax / dt);  % Numero di passi temporali
dt = Tmax / Nt;        % Ricalcolo per adattamento

% Evoluzione nel tempo con schema di Godunov
for n = 1:Nt
    u_new = u;  % Array temporaneo per aggiornamento

    for j = 2:Nx-1
        if x(j)<1
        % Flusso numerico secondo Godunov per w < 0
            F_jm12 = w * u(j);     % Flusso tra j-1 e j
            F_jp12 = w * u(j+1);   % Flusso tra j e j+1

        % Aggiornamento della soluzione
        u_new(j) = u(j) - (dt/dx) * (F_jp12 - F_jm12);
        end
    end

    % Aggiorna la soluzione
    u = u_new;

    % Plot della soluzione
    plot(x, u, 'b', 'LineWidth', 2);
    axis([-L L -0.1 2.1]);
    title(['Tempo: ', num2str(n*dt)]);
    xlabel('x'); ylabel('u(x,t)');
    drawnow;
end


