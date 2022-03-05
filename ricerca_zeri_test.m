% f = @functionToPass;
% f1 = @functionToPassDerivata;
% disp("newton");
% [x,i] = newton(f,f1,1,10^-3);
% disp(x);
% disp(i);
% [x,i] = newton(f,f1,1,10^-6);
% disp(x);
% disp(i);
% [x,i] = newton(f,f1,1,10^-9);
% disp(x);
% disp(i);
% [x,i] = newton(f,f1,1,10^-12);
% disp(x);
% disp(i);
% disp("secanti");
% [x,i] = secanti(f,1,0.99,10^-3);
% disp(x);
% disp(i);
% [x,i] = secanti(f,1,0.99,10^-6);
% disp(x);
% disp(i);
% [x,i] = secanti(f,1,0.99,10^-9);
% disp(x);
% disp(i);
% [x,i] = secanti(f,1,0.99,10^-12);
% disp(x);
% disp(i);
% disp("steffensen");
% [x,i] = steffensen(f,1,10^-3);
% disp(x);
% disp(i);
% [x,i] = steffensen(f,1,10^-6);
% disp(x);
% disp(i);
% [x,i] = steffensen(f,1,10^-9);
% disp(x);
% disp(i);
% [x,i] = steffensen(f,1,10^-12);
% disp(x);
% disp(i);
% 
f = @functionToPass2;
f1 = @functionToPass2Derivata;
disp("newton2");
[x,i] = newton(f,f1,1,10^-3);
disp(x);
disp(i);
[x,i] = newton(f,f1,1,10^-6);
disp(x);
disp(i);
[x,i] = newton(f,f1,1,10^-9);
disp(x);
disp(i);
[x,i] = newton(f,f1,1,10^-12);
disp(x);
disp(i);
disp("secanti2");
[x,i] = secanti(f,1,0.99,10^-3);
disp(x);
disp(i);
[x,i] = secanti(f,1,0.99,10^-6);
disp(x);
disp(i);
[x,i] = secanti(f,1,0.99,10^-9);
disp(x);
disp(i);
[x,i] = secanti(f,1,0.99,10^-12);
disp(x);
disp(i);
disp("steffensen2");
[x,i] = steffensen(f,1,10^-3,100000);
disp(x);
disp(i);
[x,i] = steffensen(f,1,10^-6,100000);
disp(x);
disp(i);
[x,i] = steffensen(f,1,10^-9,100000);
disp(x);
disp(i);
[x,i] = steffensen(f,1,10^-12,100000);
disp(x);
disp(i);

function y = functionToPass(x)   
    y = x - cos((pi/2)*x);
end

function y = functionToPassDerivata(x)   
    y = 1 + (pi/2)*sin((pi/2)*x);
end

function y = functionToPass2(x)   
    y = (x - cos((pi/2)*x))^3;
end

function y = functionToPass2Derivata(x)   
    y = (3*(x-cos((pi/2)*x))^2)*(1 + (pi/2)*sin((pi/2)*x));
end

function [x,i] = newton(f,f1,x0,tol,maxiter)
% [x,i] = newton(f,f1,x0,tol,maxiter)
% 
% la seguente funzione implementa il metodo di newton per 
% ricercare gli zeri di una funzione
% INPUT:
%     f - la funzione di cui si vogliono conoscere uno zero
%     f1 - la derivata prima della funzione f
%     x0 - il punto di partenza per la ricerca dello zero
%     tol - la tolleranza del metodo
%     maxiter - il numero massimo di iterazioni voluto (default 200)
% OUTPUT:
%     x - la soluzione trovata
%     i - numero di iterazioni effettuate
    if nargin == 4, maxiter = 200; end
    if maxiter < 200, maxiter = 200; end
    x_old = x0;
    for i = 0:maxiter
        x = x_old-feval(f,x_old)/feval(f1,x_old);
        if abs(x-x_old) <= tol*(1+abs(x_old)), return, end
        x_old = x;
    end
    if(i == maxiter), error("Soluzione non trovata. " + ...
            "Numero massimo di iterazioni (" + i +") raggiunto!"); end
end

function [x,i] = secanti(f,x0,x1,tol,maxiter)
% x = secanti(f,f1,x0,tol)
% 
% la funzione implementa il metdo delle secanti
% per la ricerca degli zeri di una funzione
% INPUT:
%     f - la funzione di cui si vogliono conoscere uno zero
%     x0 - il punto di partenza per la ricerca dello zero
%     x1 - il punto dopo la prima iterazione
%     tol - la tolleranza del metodo
%     maxiter - il numero massimo di iterazioni voluto (default 200)
% OUTPUT:
%     x - la soluzione trovata
%     i - numero di iterazioni effettuate
    if nargin == 4, maxiter = 200; end
    if maxiter < 200, maxiter = 200; end
    x_oldest = x0;
    x_old = x1;
    for i = 0:maxiter
        fx_old = feval(f,x_old);
        x = x_old - ((x_old-x_oldest)*fx_old)/(fx_old-feval(f,x_oldest));
        if abs(x-x_old) <= tol*(1+abs(x_old)), break; end
        appo = x_old;
        x_old = x;
        x_oldest = appo;
    end
    if(i == maxiter), print("Soluzione non trovata. " + ...
            "Numero massimo di iterazioni (" + i +") raggiunto!"); end
    return
end

function [x,i] = steffensen(f,x0,tol,maxiter)
% x = steffensen(f,f1,x0,tol,maxiter)
% 
% la funzione implementa il metdo di Steffensen
% per la ricerca degli zeri di una funzione
% INPUT:
%     f - la funzione di cui si vogliono conoscere uno zero
%     x0 - il punto di partenza per la ricerca dello zero
%     tol - la tolleranza del metodo
%     maxiter - il numero massimo di iterazioni voluto (default 200)
% OUTPUT:
%     x - la soluzione trovata
%     i - numero di iterazioni effettuate
    if nargin == 3, maxiter = 200; end
    if maxiter < 200, maxiter = 200; end
    x_old = x0;
    for i = 0:maxiter
        fx_old = feval(f,x_old);
        diff = feval(f,x_old+fx_old) - fx_old;
        if diff == 0, error("Errore: Divisone per zero! "); end
        x = x_old-(fx_old^2)/(diff);
        if abs(x-x_old) <= tol*(1+abs(x_old)), return; end
        x_old = x;
    end
    if(i == maxiter), error("Soluzione non trovata. " + ...
            "Numero massimo di iterazioni (" + i +") raggiunto!"); end
    return
end