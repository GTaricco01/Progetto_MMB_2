function dx = sirs_model(~,x,p) % t = tempo, x = [variabili], p =[parametri]
dx(1) = -p(1)*x(1)*x(2)+p(3)*x(3);
dx(2) = p(1)*x(1)*x(2) - p(2)*x(2);
dx(3) = p(2)*x(2)-p(3)*x(3);
dx = dx';
end