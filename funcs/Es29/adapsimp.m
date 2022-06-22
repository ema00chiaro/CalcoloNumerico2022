function [I4,nfeval] = adapsimp(a,b,fa,fb,tol,f)
% [I4,nfeval] = adapsimp(a,b,fa,fb,tol,f)
% 
% Calcola l'intergrale della funzione f nell'intervallo [a,b] mediante la 
% formula adattiva di Simpson
% 
% INPUT:
%     - a: estremo iniziale
%     - b: estremo finale
%     - fa: valore della funzione nel punto a
%     - fb: valore della funzione nel punto b
%     - tol: la tolleranza del metodo
%     - f: funzione integranda
% OUTPUT:
%     - I4: approssimazione dell'integrale di f nell'intervallo [a,b]
%     - nfeval : numero di valutazione funzionali effettuate

    h = (b-a)/6;
    x2 = (a+b)/2;
    f2 = feval(f,x2);
    I2 = h*(fa+4*f2+fb);
    x1 = (a+x2)/2;
    x3 = (x2+b)/2;
    f1 = feval(f,x1);
    f3 = feval(f,x3);
    I4 = I2/2 + (2*(f1+f3) - f2)*h;
    nfeval = 3;
    err = abs(I4-I2)/15;
    if err > tol
        [I4sx,nfevalsx] = adapsimp(a,x2,fa,f2,tol/2,f);
        [I4dx,nfevaldx] = adapsimp(x2,b,f2,fb,tol/2,f);
        I4 = I4sx + I4dx;
        nfeval = nfevalsx + nfevaldx;
    end
    return
end