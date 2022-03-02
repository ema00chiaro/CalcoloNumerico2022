a = [5,6,7;10,20,23;15,50,59];
b = [3;5;9];
disp(a)
a = fattLU(a);
disp(a);
x = risolviLU(a,b);
disp(x);

function a = fattLU(a)
% a = fattLU(a)
% Fattorizza la matrice in ingresso n*n in una matrice LU
% con L tr inf diag uni
% con U tr sup
% Input:
%     a - la matrice n*n
% Output:
%     a - la matrice fattorizzata LU
n= size(a);
for i = 1:n-1
    if a(i,i) == 0, error("matrice è singolare");end
    a(i+1:n, i) = a(i+1:n, i)/a(i,i); % Gauss
    a(i+1:n,i+1:n) = a(i+1:n,i+1:n) - a(i+1:n, i)*a(i, i+1:n);
end
end

function x = risolviLU(a,b)
% x = risolviLU(a,b)
% risolve un sistema la cui matrice del sistema è fattorizzata LU
% Input:
%     a - la matrice dei coefficenti
%     b - il vettore dei termini noti
% Output:
%     x - il vettore soluzione
x = b;
n = size(a);
%risolvo Ly=b
for i = 2:n
    x(i:n) = x(i:n) - a(i:n,i-1)*x(i-1);
end
%risolvo Ux=y
for i = n:-1:1
    x(i) = x(i)/a(i,i);
%     x(1:i-1) = x(1:i-1) - a(1:i-1,i)*x(i);
    x(i-1:-1:1) = x(i-1:-1:1) - a(i-1:-1:1,i)*x(i);
end
return
end