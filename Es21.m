addpath("./Utility");
a = -2;
b = 3;
n = 500;
xq = linspace(a,b,1001);
equidistanti = linspace(a,b,n+1);
yeq = functionToPass(equidistanti);
cheb = chebyshev(a,b,n);
ycheb = functionToPass(cheb);

% yqe = lagrange(equidistanti,yeq,xq);
% yqc = lagrange(cheb,ycheb,xq);

disp(interpolErrLagrange(equidistanti,yeq,xq,@functionToPass));
disp(interpolErrLagrange(cheb,ycheb,xq,@functionToPass));
disp(interpolErrNewton(equidistanti,yeq,xq,@functionToPass));
disp(interpolErrNewton(cheb,ycheb,xq,@functionToPass));

function y = functionToPass(x)
    y = 1./(2*(2*x.^2 - 2*x + 1));
    return
end

function xc = chebyshev(a,b,n)
    teta = (n:-1:0);
    teta = (teta.*2 + 1).*(pi/(2*(n+1)));
    xc = cos(teta).*((b-a)/2) + ((a+b)/2);
    return
end

function e = interpolErrLagrange(x,y,xq,f)
    e = 1 + lebesgue(x,xq);
    e = e * norm(abs(feval(f,xq)-lagrange(x,y,xq)),"inf");
end

function e = interpolErrNewton(x,y,xq,f)
    e = 1 + lebesgue(x,xq);
    e = e * norm(abs(feval(f,xq)-newton(x,y,xq)),"inf");
end

function cLebesgue = lebesgue(x,xq)
    cLebesgue = zeros(size(xq));
    for i = 1:length(x)
        cLebesgue = cLebesgue + abs(Lin(x,xq,i));
    end
    cLebesgue = norm(cLebesgue,"inf");
    return
end

function L = Lin(x,xq,i)
% L = Lin(x,xq,i)
% 
% Calcola nei punti di xq il polinomio di base Lagrange 
% di indice i definito sulle ascisse contenute in x
% 
% Input:
%     x - le ascisse sulle quali definire il polinomio
%     xq - i punti in cui si vuole calore il polinomio
%     i - indice del polinomio di base Lagrange
% Output:
%     L - il polinomio di base Lagrange calcolato nei punti xq
    n = length(x)-1;
    xi = x(i);
    x = x([1:i-1,i+1:n+1]);
    L = ones(size(xq));
    for j=1:n
        L = L.*(xq-x(j))/(xi-x(j));
    end
    return
end