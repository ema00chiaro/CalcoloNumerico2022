function [x,i] = secanti(f,x0,x1,tol,maxiter)
% x = secanti(f,f1,x0,tol)
% 
% la funzione implementa il metdo delle secanti
% per la ricerca degli zeri di una funzione
% Input:
%     f - la funzione di cui si vuole conoscere una radice
%     x0 - il punto di innesco
%     x1 - il secondo punto di innesco
%     tol - la tolleranza del metodo
%     maxiter - il numero massimo di iterazioni voluto (default 200)
% Output:
%     x - la soluzione trovata
%     i - numero di iterazioni effettuate
    if nargin == 4, maxiter = 200; end
    if feval(f,x0) == 0, i = 1; x = x0; return, end
    if feval(f,x1) == 0, i = 1; x = x1; return, end
    x_oldest = x0;
    fx_oldest = feval(f,x_oldest);
    x_old = x1;
    for i = 0:maxiter
        fx_old = feval(f,x_old);
        x = x_old - ((x_old-x_oldest)*fx_old)/(fx_old-fx_oldest);
        if abs(x-x_old) <= tol*(1+abs(x_old)), return, end
        x_oldest = x_old;
        x_old = x;
        fx_oldest = fx_old;
    end
    if(i == maxiter), error("Soluzione non trovata. " + ...
            "Numero massimo di iterazioni (" + i +") raggiunto!"); end
end

