function dx = sir_model(~,x,p) % x = [variabili], p =[parametri]
dx(1) = -p(1)*x(1)*x(2);
dx(2) = p(1)*x(1)*x(2) - p(2)*x(2);
dx(3) = p(2)*x(2);
dx = dx';
end