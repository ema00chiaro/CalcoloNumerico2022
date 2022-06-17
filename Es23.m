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

function [x,nr] = miaqr(A,b)
% [x,nr] = miaqr(A,b)
% 
% calcola la soluzione nel senso dei minimi 
% quadrati del sistema lineare Ax=b
% Input:
%     A - la matrice del sistema m x n con m >= n = rank(A)
%     b - il vettore dei termini noti di lunghezza m
% Output:
%     x - la soluzione nel senso dei minimi quadrati
%     nr - la norma del corrispondente vettore residuo
    [m,n]=size(A);
    if m < n, error("Dimensioni errate per la matrice del sistema"); end
    [row,col] = size(b);
    if col ~= 1 || row~=m, error("Il vettore dei termini noti non ha " + ...
                                 "dimensioni adeguate, errore!"); end
    for i = 1:n
        %calcolo alpha con il suo segno
        alpha = norm(A(i:m,i));
        if A(i,i) >= 0, alpha = -alpha; end
        %primo elemento di v
        v1 = A(i,i)-alpha;
        A(i+1:m,i) = A(i+1:m,i)/v1; %memorizzo il vettore v cappuccio
        A(i,i) = alpha;
        %prodotto HA -> HA = A-beta*v*vT*A
        beta = -v1/alpha;
        A(i:m,i+1:n) = A(i:m,i+1:n) - (beta*[1;A(i+1:m,i)])*([1;A(i+1:m,i)]'*A(i:m,i+1:n));
    end
    
    %calcolo g = Q'b
    g = b;
    %Hn.....H1 * b
    for i = 1:n
        v=[1;A(i+1:m,i)];
        beta = 2/(v'*v);
        g(i:m) = g(i:m) - (beta * v) * (v'*g(i:m));
    end
    
    %risolvo il sistema Rhat x = g1
    x = g(1:n);
    for i = n:-1:1
        x(i) = x(i)/A(i,i);
        x(1:i-1) = x(1:i-1) - A(1:i-1,i)*x(i);
    end

    %calcolo la norma del vettore residuo
    nr = norm(A*x-b);
    return
end
