[A,b] = linsis(10,1,1);
disp(mialdl(A,b));
[A,b] = linsis(10,10,1);
disp(mialdl(A,b));

% A = rand(5,5)*10;
% A = A'*A; %serve per renderla sdp
% disp(A);
% b = rand(5,1)*10;
% disp(b);
% disp("mat") 
% disp(A\b);
% disp("mia") 
% disp(mialdl(A,b));

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
        v = A(j,1:j-1)'.*diag(A(1:j-1,1:j-1));
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