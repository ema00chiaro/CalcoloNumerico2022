function [If,err,nfeval] = composita(fun,a,b,n,tol)
% [If,err,nfeval] = composita(fun,a,b,n,tol)
% 
% Funzione che calcola la stima err dell'errore di quadratura,
% l'approssimazione If dell'integrale della funzione fun e il
% numero di valutazioni funzionali effettuate mediante l'utilizzo
% della forma composita di Newton-Cotes base di grado n
% Input:
%     fun - la funzione integranda
%     a,b - estermi di integrazione
%     n - grado della formula composita
%     tol - tolleranza del metodo
% Output:
%     If - approssimazione dell'integrale di fun
%     err - stima dell'errore di quadratura
%     nfeval - numero di valutazioni funzionali effettuate

    if b <= a, error("a deve essere minore di b");end
    if n <= 0, error("il grado della formula" + ...
            " composita deve essere postivo");end
    nfeval = 0;
    w = weights(n);
    m = n;
    while 1
        m = 2*m;
        h = (b-a)/m;
        xi = a+(0:m)*h;
        if m == 2*n % prima volta
            fi = fun(xi);
            nfeval = nfeval + m+1;
        else
            fi = repelem(fi,2);
            fi = fi(1:m+1);
            fi(2:2:m) = fun(xi(2:2:m));
            nfeval = nfeval + m/2;
        end
        if mod(n,2) == 0
            mu = 2;
        else
            mu = 1;
        end
        
        Ieven = 0;
        i = 1;
        k = 1;
        while i <= m+1
            Ieven = Ieven + w(k)*fi(i);
            i = i+2;
            k = k+1;
            if k == n+2
                k = 1;
                i = i-2;
            end
        end
        Ieven = Ieven * (b-a)/(m/2);
        
        If = 0;
        i = 1;
        k = 1;
        while i <= m+1
            If = If + w(k)*fi(i);
            i = i+1;
            k = k+1;
            if k == n+2
                k = 1;
                i = i-1;
            end
        end
        If = If * (b-a)/m;
    
        err = abs(Ieven-If)/(2^(n+mu)-1);
        if err <= tol,break;end
    end
end