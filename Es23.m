n = 10000;
xi = (0:n)./n;
yi = functionToPassPerturb(xi);
a = 0;
b = 1;
m = 15;
% p = polyfit(xi,yi,m);
% xq = linspace(a,b,1000);
approxErrorGraph(a,b,m,n,@functionToPass,@functionToPassPerturb);

function approxErrorGraph(a,b,m,n,f,fpert)
% approxErrorGraph(a,b,m,n,f,fpert)
% 
% Funzione che grafica in formato semilogy l'errore di 
% approssimazione ||f - p||, stimato come il massimo errore sui punti di xi
% equidistanti fra loro, all'interno di un intervallo [a,b] 
% per i vari polinomi interpolanti fpert,
% una perturbazione di f, di grado da 1 a m
% 
% Input:
%     a,b - gli estremi dell'intervallo
%     m - il grado massimo del polinomio interpolante fpert
%     n - il quantitativo delle ascisse di interpolazione
%     f - la funzione non perturbata
%     fpert - la funzione perturbata
% Output:
%     Il grafico semilogy di ||f-p|| per i valori da 1 a m
    
    if b <= a, error("a deve essere minore di b");end
    if m <= 0, error("il grado del polinomio deve " + ...
            "essere positivo");end
    if n <= 0, error("numero delle ascisse di interpolazione" + ...
            " deve essere positivo");end
 
    xi = (0:n)./n;
    yi = feval(fpert,xi);
    e = zeros(1,m);
    for i = 1:m
        p = polyfit(xi,yi,i);
        e(i) = norm(f(xi)-polyval(p,xi),"inf");
    end
    x_axis = (1:m);
    semilogy(x_axis,e);
    xlabel("m");
    ylabel("e = || f-p ||");
end


% plot(xq,polyval(p,xq));
% hold on;
% plot(xq,functionToPass(xq));
% hold off;

function y = functionToPass(x)
    y = sin(pi*x.^2);
end

function y = functionToPassPerturb(x)
y = functionToPass(x) + 1e-1 * rand(size(x));
end
