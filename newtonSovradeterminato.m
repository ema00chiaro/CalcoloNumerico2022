format long;
fun=@(x)[(x(1)^2+1)*(x(2)-2); exp(x(1)-1)+exp(x(2)-2)-2];
Df=@(x)[2*x(1)*(x(2)-2) (x(1)^2+1);exp(x(1)-1) exp(x(2)-2)];

fun2=@(x)[x(1)-(x(2)*x(3));exp(x(1)+x(2)+x(3)-3)-x(2);x(1)+x(2)+2*x(3)-4];
Df2=@(x)[1 -x(3) -x(2);exp(x(1)+x(2)+x(3)-3) exp(x(1)+x(2)+x(3)-3)-1 exp(x(1)+x(2)+x(3)-3);1 1 2];
[x,nit]=newton(fun,Df,[0,0]',1e-3,200);
disp("iter: " + nit);
disp(x);
[x,nit]=newton(fun,Df,[0,0]',1e-8,200);
disp("iter: " + nit);
disp(x);
[x,nit]=newton(fun,Df,[0,0]',1e-13,200);
disp("iter: " + nit);
disp(x);
[x2,nit2]=newton(fun2,Df2,[0,0,0]',1e-3,200);
disp("iter: " + nit2);
disp(x2);
[x2,nit2]=newton(fun2,Df2,[0,0,0]',1e-8,200);
disp("iter: " + nit2);
disp(x2);
[x2,nit2]=newton(fun2,Df2,[0,0,0]',1e-13,200);
disp("iter: " + nit2);
disp(x2);
function [x,nit] = newton(fun,jacobian,x0,tol,maxit)
% [x,nit] = newton(fun,jacobian,x0,tol,maxit)
% 
% La funzione risolve un sistema di equazioni non lineari
% mediante il metodo di newton.
% Input:
%     fun - sistema di equazioni da risolvere
%     jacobian - matrice jacobiana del sistema
%     x0 - punto di innesco per la ricerca della soluzione
%     tol - tolleranza del metodo
%     maxit - numero di iterazioni massime desiderate
% Output:
%     x - soluzione del sistema non lineare
%     nit - numero di iterazioni effettuate
    xold=x0;
    deltax=mialu(jacobian(xold),-fun(xold));
    xnew=xold+deltax;
    nit=1;
    while norm((xnew-xold)./(1+abs(xold)))>tol && nit<=maxit 
        xold=xnew;
        deltax=mialu(jacobian(xold),-fun(xold));
        xnew=xold+deltax;
        nit=nit+1;
    end
    x=xnew;
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