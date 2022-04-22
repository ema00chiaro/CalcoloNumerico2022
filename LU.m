% [A,b] = linsis(10,1);
% disp(mialu(A,b));
% [A,b] = linsis(10,10);
% disp(cond(A));
% disp(mialu(A,b));

A = rand(5,5)*100;
disp(A);
b = rand(5,1)*10;
disp(b);
disp("mat") 
disp(A\b);
disp("mia") 
disp(mialu(A,b));

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
    if col ~= 1 || row ~= n
        error("Il vettore dei termini noti " + ...
            "non ha dimensioni adeguate, errore!");
    end
    for i = 1:n-1
        if A(i,i) == 0, error("Matrice singolare"); end
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