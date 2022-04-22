format long;
fun=@(x)[(x(1)^2+1)*(x(2)-2); exp(x(1)-1)+exp(x(2)-2)-2];
Df=@(x)[2*x(1)*(x(2)-2) (x(1)^2+1);exp(x(1)-1) exp(x(2)-2)];

fun2=@(x)[x(1)-(x(2)*x(3));exp(x(1)+x(2)+x(3)-3)-x(2);x(1)+x(2)+2*x(3)-4];
Df2=@(x)[1 -x(3) -x(2);exp(x(1)+x(2)+x(3)-3) exp(x(1)+x(2)+x(3)-3)-1 exp(x(1)+x(2)+x(3)-3);1 1 2];
[x,nit]=newton(fun,Df,[0,0]',1e-3,200);
disp(x + ":" + nit);
[x,nit]=newton(fun,Df,[0,0]',1e-8,200);
disp(x + ":" + nit);
[x,nit]=newton(fun,Df,[0,0]',1e-13,200);
disp(x + ":" + nit);
[x2,nit2]=newton(fun2,Df2,[0,0,0]',1e-3,200);
disp(x2 + ":" + nit2);
[x2,nit2]=newton(fun2,Df2,[0,0,0]',1e-8,200);
disp(x2 + ":" + nit2);
[x2,nit2]=newton(fun2,Df2,[0,0,0]',1e-13,200);
disp(x2 + ":" + nit2);
function [x,nit] = newton(fun,jacobian,x0,tol,maxit)
% [x,nit] = newton(fun,jacobian,x0,tol,maxit)
%
% La funzione risolve un sistema di equazioni non lineari.
% fun - sistema di equazioni da risolvere
% jacobian - matrice jacobiana del sistema dato in ingresso
% x0 - punto di inizio per la ricerca della soluzione
% tol - tolleranza del metodo
% maxit - numero di iterazione massime del metodo
xold=x0;
deltax=palu(jacobian(xold),-fun(xold));
xnew=xold+deltax;
nit=1;
while norm((xnew-xold)./(1+abs(xold)))>tol && nit<=maxit 
    xold=xnew;
    deltax=palu(jacobian(xold),-fun(xold));
    xnew=xold+deltax;
    nit=nit+1;
end
x=xnew;
end



function [a,pos] = fattPALU(a)
    n = size(a);
    pos = (1:n)';
    for i = 1:n-1
        %prendo il massimo della colonna in valore assoluto
        [mi,ki] = max(abs(a(i:n,i)));
        if mi == 0, error("matrice singolare");end
        ki = ki+i-1;
        %scambio le righe
        pos([i,ki])= pos([ki,i]);
        a([i,ki],:)=a([ki,i],:);
        %metto la colonna da i+1 tutta a zero figuaratamente, metto in quella colonna il vettore di gauss
        a(i+1:n,i) = a(i+1:n,i)/a(i,i);
        %Ai+1 = Ai-gi*eiT*Ai
        a(i+1:n,i+1:n) = a(i+1:n,i+1:n) - a(i+1:n,i)*a(i,i+1:n);
    end
    return
end

function x = solvePALU(a,b,pos)
    x = b(pos);
    n = size(a);
    %Ly=b(pos)
    for i = 2:n
        x(i:n) = x(i:n) - a(i:n,i-1)*x(i-1);
    end
    %Ux=y
    for i = n:-1:1
        x(i) = x(i)/a(i,i);
        x(1:i-1) = x(1:i-1) - a(1:i-1,i)*x(i);
    end
    return
end

function x = palu(a,b)
    n = size(a);
    pos = (1:n)';
    for i = 1:n-1
        %prendo il massimo della colonna in valore assoluto
        [mi,ki] = max(abs(a(i:n,i)));
        if mi == 0, error("matrice singolare");end
        ki = ki+i-1;
        %scambio le righe
        pos([i,ki])= pos([ki,i]);
        a([i,ki],:)=a([ki,i],:);
        %metto la colonna da i+1 tutta a zero figuaratamente, metto in quella colonna il vettore di gauss
        a(i+1:n,i) = a(i+1:n,i)/a(i,i);
        %Ai+1 = Ai-gi*eiT*Ai
        a(i+1:n,i+1:n) = a(i+1:n,i+1:n) - a(i+1:n,i)*a(i,i+1:n);
    end

    x = b(pos);
    %Ly=b(pos)
    for i = 2:n
        x(i:n) = x(i:n) - a(i:n,i-1)*x(i-1);
    end
    %Ux=y
    for i = n:-1:1
        x(i) = x(i)/a(i,i);
        x(1:i-1) = x(1:i-1) - a(1:i-1,i)*x(i);
    end
end