function [u] = Godunov(L,Nx,Tmax,CFL,w,u_0)

%u_0 deve essere un vettore che mi dà le condizoni iniziali del ciclo che
%voglo mettere in atto.

% Discretizzazione spaziale
dx = 2*L / Nx;
x = linspace(-L, L, Nx);

% Calcolo del passo temporale (Condizione CFL corretta)
dt = CFL * (dx /abs(w));
Nt = ceil(Tmax / dt);  % Numero di passi temporali
dt = Tmax / Nt;        % Ricalcolo per adattamento

% Evoluzione nel tempo con schema di Godunov
u=u_0;
for n = 1:Nt
    u_new = u;  % Array temporaneo per aggiornamento

    for j = 2:Nx-1
        if abs(x(j))<1
        % Flusso numerico secondo Godunov per w < 0
            F_jm12 = w * u(j);     % Flusso tra j-1 e j
            F_jp12 = w * u(j+1);   % Flusso tra j e j+1

        % Aggiornamento della soluzione
        u_new(j) = u(j) - (dt/dx) * (F_jp12 - F_jm12);
        end
    end

    % Aggiorna la soluzione
    u = u_new;

end
