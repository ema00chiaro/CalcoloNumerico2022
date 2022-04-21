function [x,nit] = newton(fun,jacobian,x0,tol,maxit)
%% La funzione risolve un sistema di equazioni non lineari.
%% fun - sistema di equazioni da risolvere
%% jacobian - matrice jacobiana del sistema dato in ingresso
%% x0 - punto di inizio per la ricerca della soluzione
%% tol - tolleranza del metodo
%% maxit - numero di iterazione massime del metodo
%%
xold=x0;
xnew=x0-(jacobian(x0))\(fun(x0));
nit=1;
while norm((xnew-xold)./(1+abs(xold)))>tol && nit<=maxit 
    xold=xnew;
    xnew=xold-(jacobian(xold))\fun(xold);
    nit=nit+1;
end
x=xnew;
display(x);
end

