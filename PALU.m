l = [1,0,0,0;1/2,1,0,0;3/4,1/2,1,0;0,1/2,1/2,1];
u=[4,0,1,1;0,2,7/2,1/2;0,0,1/2,0;0,0,0,-1/4];
a=[4,0,1,1;3,1,3,1;0,1,2,0;2,2,4,1];
p=[1,0,0,0;0,0,0,1;0,1,0,0;0,0,1,0];
b=[3;5;9;6];
[a,pos] = fattPALU(a);
disp(a);
x = solvePALU(a,b,pos);
disp(x);


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