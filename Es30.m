tols = [10^-2;10^-3;10^-4;10^-5;10^-6];
tol = 10^-2;
ris2 = zeros(1,5);
I2 = adaptrap(0,1,f(0),f(1),tol,@f);
ris2(1) = I2;
tol = 10^-3;
I2 = adaptrap(0,1,f(0),f(1),tol,@f);
ris2(2) = I2;
tol = 10^-4;
I2 = adaptrap(0,1,f(0),f(1),tol,@f);
ris2(3) = I2;
tol = 10^-5;
I2 = adaptrap(0,1,f(0),f(1),tol,@f);
ris2(4) = I2;
tol = 10^-6;
I2 = adaptrap(0,1,f(0),f(1),tol,@f);
ris2(5) = I2;

trap = ris2';

ris4 = zeros(1,5);
tol = 10^-2;
I4 = adapsimp(0,1,f(0),f(1),tol,@f);
ris4(1) = I4;
tol = 10^-3;
I4 = adapsimp(0,1,f(0),f(1),tol,@f);
ris4(2) = I4;
tol = 10^-4;
I4 = adapsimp(0,1,f(0),f(1),tol,@f);
ris4(3) = I4;
tol = 10^-5;
I4 = adapsimp(0,1,f(0),f(1),tol,@f);
ris4(4) = I4;
tol = 10^-6;
I4 = adapsimp(0,1,f(0),f(1),tol,@f);
ris4(5) = I4;

simp = ris4';
disp(table(tols,trap,simp));

function y = f(x)   
    y = sin(1/(0.01+x));
end

function I2 = adaptrap(a,b,fa,fb,tol,f)
% I2 = adaptrap(a,b,fa,fb,tol,f)
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

    h = b-a;
    x1 = (a+b)/2;
    f1 = feval(f,x1);
    I1 = (h/2)*(fa+fb);
    I2 = (I1 + h*f1)/2;
    err = abs(I2-I1)/3;
    if err > tol
        I2 = adaptrap(a,x1,fa,f1,tol/2,f) + adaptrap(x1,b,f1,fb,tol/2,f);
    end
    return
end

function I4 = adapsimp(a,b,fa,fb,tol,f)
% I4 = adapsimp(a,b,fa,fb,tol,f)
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

    h = (b-a)/6;
    x2 = (a+b)/2;
    f2 = feval(f,x2);
    I2 = h*(fa+4*f2+fb);
    x1 = (a+x2)/2;
    x3 = (x2+b)/2;
    f1 = feval(f,x1);
    f3 = feval(f,x3);
    I4 = I2/2 + (2*(f1+f3) - f2)*h;
    err = abs(I4-I2)/15;
    if err > tol
        I4 = adapsimp(a,x2,fa,f2,tol/2,f) + adapsimp(x2,b,f2,fb,tol/2,f);
    end
    return
end