function [x,nit] = newton(fun,jacobian,x0,tol,maxit)
% [x,nit] = newton(fun,jacobian,x0,tol,maxit)
% 
% La funzione risolve un sistema di equazioni non lineari
% mediante il metodo di newton.
% Input:
%     fun - sistema di equazioni da risolvere
%     jacobian - matrice jacobiana del sistema
%     x0 - punto di innesco per la ricerca della soluzione
%     tol - tolleranza del metodo
%     maxit - numero di iterazioni massime desiderate
% Output:
%     x - soluzione del sistema non lineare
%     nit - numero di iterazioni effettuate
    xold=x0;
    deltax=mialu(jacobian(xold),-fun(xold));
    xnew=xold+deltax;
    nit=1;
    while norm((xnew-xold)./(1+abs(xold)))>tol && nit<=maxit 
        xold=xnew;
        deltax=mialu(jacobian(xold),-fun(xold));
        xnew=xold+deltax;
        nit=nit+1;
    end
    x=xnew;
end