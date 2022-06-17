% a = [1,1,1;2,10,5;1,1,1;6,7,4];
% b=[3;5;9;6];
% [x,nr] = miaqr(a,b);
% disp(x);
% disp(nr);
A = [ 1 3 2; 3 5 4; 5 7 6; 3 6 4; 1 4 2 ];
b = [ 15 28 41 33 22 ]';
D = diag(1:5);
disp(A\b);
[x,nr] = miaqr(A,b);
disp(x);
disp("err: " + abs(A\b-x)./abs(A\b));
% disp(nr);
disp(D*A\D*b);
[x,nr] = miaqr(D*A,D*b);
disp(x);
disp("err: " + abs((D*A)\(D*b)-x)./abs((D*A)\(D*b)));
% disp(nr);


function [x,nr] = miaqr(A,b)
% [x,nr] = miaqr(A,b)
% 
% calcola la soluzione nel senso dei minimi 
% quadrati del sistema lineare Ax=b
% Input:
%     A - la matrice del sistema m x n con m >= n = rank(A)
%     b - il vettore dei termini noti di lunghezza m
% Output:
%     x - la soluzione nel senso dei minimi quadrati
%     nr - la norma del corrispondente vettore residuo
    [m,n]=size(A);
    if m < n, error("Dimensioni errate per la matrice del sistema"); end
    [row,col] = size(b);
    if col ~= 1 || row~=m, error("Il vettore dei " + ...
            "termini noti non ha dimensioni adeguate, errore!"); end
    for i = 1:n
        %calcolo alpha con il suo segno
        alpha = norm(A(i:m,i));
        if A(i,i) >= 0, alpha = -alpha; end
        %primo elemento di v
        v1 = A(i,i)-alpha;
        A(i+1:m,i) = A(i+1:m,i)/v1; %memorizzo il vettore v cappuccio
        A(i,i) = alpha;
        %prodotto HA -> HA = A-beta*v*vT*A
        beta = -v1/alpha;
        A(i:m,i+1:n) = A(i:m,i+1:n) - (beta*[1;A(i+1:m,i)])*([1;A(i+1:m,i)]'*A(i:m,i+1:n));
    end
    
    %calcolo g = Q'b
    g = b;
    %Hn.....H1 * b
    for i = 1:n
        v=[1;A(i+1:m,i)];
        beta = 2/(v'*v);
        g(i:m) = g(i:m) - (beta * v) * (v'*g(i:m));
    end
    
    %risolvo il sistema Rhat x = g1
    x = g(1:n);
    for i = n:-1:1
        x(i) = x(i)/A(i,i);
        x(1:i-1) = x(1:i-1) - A(1:i-1,i)*x(i);
    end

    %calcolo la norma del vettore residuo
    nr = norm(A*x-b);
    return
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


