function dx = sir_model(~,x,p,class,M) % x = [variabili], p =[beta,gamma]
dx(1) = -p(1)*x(1)*sum((M(class,:)*x(2)));
dx(2) = p(1)*x(1)*sum((M(class,:)*x(2))) - p(2)*x(2);
dx(3) = p(2)*x(2);
dx = dx';
end