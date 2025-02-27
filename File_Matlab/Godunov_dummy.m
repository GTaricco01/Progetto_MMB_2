clearvars; clc; close all;

% Parametri
L    = 10;                               % Lunghezza del dominio
Nu   = 1000;                            % Numero di punti spaziali
Tmax = 50;                            % Tempo massimo di simulazione
CFL  = 0.9;                            % Numero di Courant-Friedrichs-Lewy
% w    = @(u) -.5*(abs(u)<=1);                             % Velocità di drift (negativa)
w    = @(u) -.5

% Discretizzazione spaziale
du = 2*L / Nu;
uu = linspace(-L, L, Nu);

 
% Condizione iniziale: onda quadrata
f0 = @(u) 2*(abs(u)<1)+ 10*(u<=-1); 
f  = zeros(1, Nu);

for i = 1:Nu-1
    f(i) = (f0(uu(i)) + f0(uu(i+1))) / 2; % Media di cella per inizializzazione
end

% Calcolo del passo temporale (Condizione CFL corretta)
dt = CFL * (du /max(abs(w(uu))));
Nt = ceil(Tmax / dt);  % Numero di passi temporali
dt = Tmax / Nt;        % Ricalcolo per adattamento

% Evoluzione nel tempo con schema di Godunov
for n = 1:Nt
    f_new = f;  % Array temporaneo per aggiornamento
    for j = 2:Nu-1
        if (uu(j))<1
        % Flusso numerico secondo Godunov per w < 0
            g_jm12 = w(uu(j)) * f(j);     % Flusso tra j-1 e j
            g_jp12 = w(uu(j)) * f(j+1);   % Flusso tra j e j+1
        % Aggiornamento della soluzione
        f_new(j) = f(j) - (dt/du) * (g_jp12 - g_jm12);
        end
    end
    % Aggiorna la soluzione
    f = f_new;
    % Plot della soluzione
    plot(uu, f, 'b', 'LineWidth', 2);
    axis([-L L -0.01 10]);
    title(['Tempo: ', num2str(n*dt)]);
    xlabel('u'); ylabel('f(u,t)');
    drawnow;
end


