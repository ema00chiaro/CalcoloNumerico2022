function x = mialdlt(A,b)
% x = mialdlt(A,b)
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
    if col ~= 1 || row ~= n
        error("Il vettore dei termini noti non ha " + ...
            "dimensioni adeguate, errore!");
    end
    
    if A(1,1) <= 0, error("Matrice non sdp");end
    A(2:n,1) = A(2:n,1)/A(1,1);
    for j = 2:n
        v = A(j,1:j-1)'.*diag(A(1:j-1,1:j-1));
        A(j,j) = A(j,j) - A(j,1:j-1)*v;
        if A(j,j) <= 0, error("Trovato un elemento negativo" + ...
                " sulla diagonale!");end
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