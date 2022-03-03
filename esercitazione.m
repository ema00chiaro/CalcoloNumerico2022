
% x = input("numero di cui trovare la radice: ");
% x = radice(x);
% disp("valore : " + x);
% testRadice();

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
% f = @functionToPass2;
% f1 = @functionToPass2Derivata;
% disp("newton2");
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
% disp("secanti2");
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
% disp("steffensen2");
% [x,i] = steffensen(f,1,10^-3,100000);
% disp(x);
% disp(i);
% [x,i] = steffensen(f,1,10^-6,100000);
% disp(x);
% disp(i);
% [x,i] = steffensen(f,1,10^-9,100000);
% disp(x);
% disp(i);
% [x,i] = steffensen(f,1,10^-12,100000);
% disp(x);
% disp(i);

[A,b] = linsis(10,1);
disp(mialu(A,b));
[A,b] = linsis(10,10);
disp(mialu(A,b));
[A,b] = linsis(10,1,1);
disp(mialdl(A,b));
[A,b] = linsis(10,10,1);
disp(mialdl(A,b));


disp("ciao");
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
        x = x_old-(fx_old^2)/(feval(f,x_old+fx_old)-fx_old);
        if abs(x-x_old) <= tol*(1+abs(x_old)), return; end
        x_old = x;
    end
    if(i == maxiter), error("Soluzione non trovata. " + ...
            "Numero massimo di iterazioni (" + i +") raggiunto!"); end
end

function x = radice(x)
% x = radice(x)
% 
% la funziona calcola la radice di un numero mediante le operazioni elementari
% INPUT:
%     x - il numero di cui si vuole trovare la radice
% OUTPUT:
%     x - la radice del numero
    if x < 0
        error("errore: numero negativo, impossibile fare la radice!"); 
    end
    if (x == 0 || x == 1), return, end
    sq = x;
    sq_old = 0;
    i = 0;
    while (abs(sq-sq_old) > eps)
        i = i+1;
        sq_old = sq;
        sq = (1/2)*(sq_old+(x/sq_old));
    end
    x = sq;
    return
end


function testRadice()
    inf = 1e-10;
    sup = 1e10;
    x = inf;
    err = 0;
    i = 0;
    while x <= sup
        i = i+1;
        mia = radice(x);
        vera = sqrt(x);
        disp("ITERAZIONE " + i);
        if (abs(vera - mia) > eps)
            err = err + 1; 
            disp("vera: "+ vera);
            disp("mia: "+ mia);
            disp ("diffe: " + abs(vera - mia));
            disp("QUI QUALCOSA NON VA!");
        end
        x = x*10;
    end
    disp("errori totali: " + err)
end

function x = mialu(A,b)
% x = mialu(A,b)
% 
% La funzione calcola la soluzione del sistema lineare Ax=b 
% mediante il metodo della fattorizzazione LU con pivoting parziale
% Input:
%     A - la matrice del sistema (n*n)
%     b - il vettore dei termini noti (n)
% Output:
%     x - la soluzione trovata
    [row,col]= size(A);
    if row ~= col, error("Matrice non quadrata, errore!"); end
    n = row;
    [row,col]= size(b);
    if col ~= 1 || row ~= n, error("Il vettore dei termini noti non ha " + ...
                            "dimensioni adeguate, errore!"); end
    for i = 1:n-1
        if A(i,i) == 0, disp("indice " + i + " valore " + A(i,i)); error("Matrice singolare");end
        A(i+1:n, i) = A(i+1:n, i)/A(i,i); % Gauss
        A(i+1:n,i+1:n) = A(i+1:n,i+1:n) - A(i+1:n, i)*A(i, i+1:n);
    end
    x = b;
    %Ly=b
    for i = 2:n
        x(i:n) = x(i:n) - A(i:n,i-1)*x(i-1);
    end
    %Ux=y
    for i = n:-1:1
        x(i) = x(i)/A(i,i);
        x(1:i-1) = x(1:i-1) - A(1:i-1,i)*x(i);
    end
    return
end

function x = mialdl(A,b)
% x = mialu(A,b)
% 
% La funzione calcola la soluzione del sistema lineare Ax=b 
% mediante il metodo della fattorizzazione LDL'
% Input:
%     A - la matrice del sistema (n*n)
%     b - il vettore dei termini noti (n)
% Output:
%     x - la soluzione trovata
    [row,col]= size(A);
    if row ~= col, error("Matrice non quadrata, errore!"); end
    n = row;
    [row,col]= size(b);
    if col ~= 1 || row ~= n, error("Il vettore dei termini noti non ha " + ...
                            "dimensioni adeguate, errore!"); end
    
    if A(1,1) <= 0, error("Matrice non sdp");end
    A(2:n,1) = A(2:n,1)/A(1,1);
    for j = 2:n
        %FIXME da rivedere questa cosa qua
        v = diag(diag(A(1:j-1,1:j-1)))*A(j,1:j-1)';
        A(j,j) = A(j,j) - A(j,1:j-1)*v;
        if A(j,j) <= 0, error("Trovato un elemento negativo sulla diagonale!");end
        A(j+1:n,j) = (A(j+1:n,j) - A(j+1:n,1:j-1)*v)/A(j,j);
    end
    
    x = b;
    %Ly2=b
    for i = 2:n
        x(i:n) = x(i:n) - A(i:n,i-1)*x(i-1);
    end
    %Dy1=y2
    x = x./diag(A);
    %L'x=y1
    for i=n-1:-1:1
        x(1:i) = x(1:i)-A(i+1,1:i)'*x(i+1);
    end
end

function [A,b] = linsis(n,k,simme)
% [A,b] = linsis(n,k,simme)
% 
% Crea una matrice A nxn ed un termine noto b,
% in modo che la soluzione del sistema lineare
% A*x=b sia x = [1,2,...,n]'.
% INPUT:
%     k - e' un parametro ausiliario.
%     simme - se specificato, crea una matrice simmetrica e definita positiva.
% OUTPUT:
%     A - La matrice nxn ritornata
%     b - Il vettore dei termini noti ritornato
    sigma = 10^(-2*(1-k))/n;
    rng(0);
    [q1,~] = qr(rand(n));
    if nargin==3
        q2 = q1';
    else
        [q2,~] = qr(rand(n));
    end
    A = q1*diag([sigma 2/n:1/n:1])*q2;
    x = (1:n)';
    b = A*x;
    return
end




























