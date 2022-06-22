function [I2,nfeval] = adaptrap(a,b,fa,fb,tol,f)
% [I2,nfeval] = adaptrap(a,b,fa,fb,tol,f)
% 
% Calcola l'intergrale della funzione f nell'intervallo [a,b] mediante la 
% formula adattiva dei trapezi
% 
% INPUT:
%     - a: estremo iniziale
%     - b: estremo finale
%     - fa: valore della funzione nel punto a
%     - fb: valore della funzione nel punto b
%     - tol: la tolleranza del metodo
%     - f: funzione integranda
% OUTPUT:
%     - I2: approssimazione dell'integrale di f nell'intervallo [a,b]
%     nfeval - numero di valutazioni funzionali effettuate
    
    h = b-a;
    x1 = (a+b)/2;
    f1 = feval(f,x1);
    nfeval = 1;
    I1 = (h/2)*(fa+fb);
    I2 = (I1 + h*f1)/2;
    err = abs(I2-I1)/3;
    if err > tol
        [I2sx,nfevalsx] = adaptrap(a,x1,fa,f1,tol/2,f);
        [I2dx,nfevaldx] = adaptrap(x1,b,f1,fb,tol/2,f);
        I2 = I2sx + I2dx;
        nfeval = nfevalsx + nfevaldx;
    end
    return
end