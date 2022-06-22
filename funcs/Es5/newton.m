function [x,i] = newton(f,f1,x0,tol,maxiter)
% [x,i] = newton(f,f1,x0,tol,maxiter)
% 
% la seguente funzione implementa il metodo di newton per 
% ricercare gli zeri di una funzione
% Input:
%     f - la funzione di cui si vuole conoscere una radice
%     f1 - la derivata prima della funzione f
%     x0 - il punto di innesco
%     tol - la tolleranza del metodo
%     maxiter - il numero massimo di iterazioni voluto (default 200)
% Output:
%     x - la soluzione trovata
%     i - numero di iterazioni effettuate
    if nargin == 4, maxiter = 200; end
    if feval(f,x0) == 0, i = 1; x = x0; return, end
    x_old = x0;
    for i = 0:maxiter
        x = x_old-feval(f,x_old)/feval(f1,x_old);
        if abs(x-x_old) <= tol*(1+abs(x_old)), return, end
        x_old = x;
    end
    if(i == maxiter), error("Soluzione non trovata. " + ...
            "Numero massimo di iterazioni (" + i +") raggiunto!"); end
end