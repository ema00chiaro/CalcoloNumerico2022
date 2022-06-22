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
    pos = (1:n)';
    for i = 1:n-1
        %prendo il massimo della colonna in valore assoluto
        [mi,ki] = max(abs(A(i:n,i)));
        if mi == 0, error("matrice singolare");end
        ki = ki+i-1;
        %scambio le righe
        pos([i,ki])= pos([ki,i]);
        A([i,ki],:)=A([ki,i],:);
        %utilizzo questa colonna per inserire il vettore di Gauss
        A(i+1:n,i) = A(i+1:n,i)/A(i,i);
        %Ai+1 = Ai-gi*eiT*Ai
        A(i+1:n,i+1:n) = A(i+1:n,i+1:n) - A(i+1:n,i)*A(i,i+1:n);
    end

    x = b(pos);
    %Ly=b(pos)
    for i = 2:n
        x(i:n) = x(i:n) - A(i:n,i-1)*x(i-1);
    end
    %Ux=y
    for i = n:-1:1
        x(i) = x(i)/A(i,i);
        x(1:i-1) = x(1:i-1) - A(1:i-1,i)*x(i);
    end
end