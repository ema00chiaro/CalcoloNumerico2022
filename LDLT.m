a = [6,24,30;24,97,126;30,126,217];
b= [3;5;9];

a = fattLDLT(a);
disp(a);
x = solveLDLT(a,b);
disp(x);

function a = fattLDLT(a)
    n = size(a);
    if a(1,1) <= 0, error("matrice non sdp");end
    a(2:n,1) = a(2:n,1)/a(1,1);
    for j = 2:n
        v = diag(diag(a(1:j-1,1:j-1)))*a(j,1:j-1)';
        a(j,j) = a(j,j) - a(j,1:j-1)*v;
        if a(j,j) <= 0, error("non ha diagonale positiva");end
        a(j+1:n,j) = (a(j+1:n,j) - a(j+1:n,1:j-1)*v)/a(j,j);
    end
end

function x = solveLDLT(a,b)
    x = b;
    n = size(a);
    %Ly2=b
    for i = 2:n
        x(i:n) = x(i:n) - a(i:n,i-1)*x(i-1);
    end
    %Dy1=y2
    x = x./diag(a);
    %LTx=y1
    for i=n-1:-1:1
        x(1:i) = x(1:i)-a(i+1,1:i)'*x(i+1);
    end
end