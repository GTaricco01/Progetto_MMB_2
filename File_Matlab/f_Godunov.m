function f = f_Godunov(L,Nu,Tmax,CFL,w,f_0)
% f_0 deve essere un vettore che mi dà le condizoni iniziali del ciclo che
% voglo mettere in atto.

% Discretizzazione spaziale
du = 2*L / Nu;
uu = linspace(-L, L, Nu);
% Calcolo del passo temporale (Condizione CFL corretta)
dt = CFL * (du /abs(w));
Nt = ceil(Tmax / dt);  % Numero di passi temporali
dt = Tmax / Nt;        % Ricalcolo per adattamento
% Evoluzione nel tempo con schema di Godunov
f = f_0;
for n = 1:Nt
    f_new = f;  % Array temporaneo per aggiornamento
    for j = 2:Nu-1
        if abs(uu(j))<1
        % Flusso numerico secondo Godunov per w < 0
            g_jm12 = w(uu(j)) * f(j);     % Flusso tra j-1 e j
            g_jp12 = w(uu(j)) * f(j+1);   % Flusso tra j e j+1
        % Aggiornamento della soluzione
        f_new(j) = f(j) - (dt/du) * (g_jp12 - g_jm12);
        end
    end
    % Aggiorna la soluzione
    f = f_new;
end
