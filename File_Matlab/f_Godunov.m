function f = f_Godunov(L,Nx,Tmax,CFL,w,f_0)

% f_0 deve essere un vettore che mi dà le condizoni iniziali del ciclo che
% voglo mettere in atto.

% Discretizzazione spaziale
dx = 2*L / Nx;
x = linspace(-L, L, Nx);

% Calcolo del passo temporale (Condizione CFL corretta)
dt = CFL * (dx /abs(w));
Nt = ceil(Tmax / dt);  % Numero di passi temporali
dt = Tmax / Nt;        % Ricalcolo per adattamento

% Evoluzione nel tempo con schema di Godunov
f = f_0;
for n = 1:Nt
    f_new = f;  % Array temporaneo per aggiornamento

    for j = 2:Nx-1
        if abs(x(j))<1
        % Flusso numerico secondo Godunov per w < 0
            g_jm12 = w * f(j);     % Flusso tra j-1 e j
            g_jp12 = w * f(j+1);   % Flusso tra j e j+1

        % Aggiornamento della soluzione
        f_new(j) = f(j) - (dt/dx) * (g_jp12 - g_jm12);
        end
    end

    % Aggiorna la soluzione
    f = f_new;

end
