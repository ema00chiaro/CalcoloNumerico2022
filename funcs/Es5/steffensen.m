function [x,i] = steffensen(f,x0,tol,maxiter)
% x = steffensen(f,f1,x0,tol,maxiter)
% 
% la funzione implementa il metodo di Steffensen
% per la ricerca degli zeri di una funzione
% Input:
%     f - la funzione di cui si vuole conoscere una radice
%     x0 - il punto di innesco
%     tol - la tolleranza del metodo
%     maxiter - il numero massimo di iterazioni voluto (default 200)
% Output:
%     x - la soluzione trovata
%     i - numero di iterazioni effettuate

    if nargin == 3, maxiter = 200; end
    if feval(f,x0) == 0, i = 1; x = x0; return, end
    x_old = x0;
    for i = 0:maxiter
        fx_old = feval(f,x_old);
        diff = feval(f,x_old+fx_old) - fx_old;
        if diff == 0, error("Errore: Divisone per zero! -Iterazione: " + (i+1)); end
        x = x_old-(fx_old^2)/diff;
        if abs(x-x_old) <= tol*(1+abs(x_old)), return, end
        x_old = x;
    end
    if(i == maxiter), error("Soluzione non trovata. " + ...
            "Numero massimo di iterazioni (" + i +") raggiunto!"); end
end

