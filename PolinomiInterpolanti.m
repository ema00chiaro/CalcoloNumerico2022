% x = [1,1,2,2,3,3];
% y = [11,10,23,22,34,33];
% disp(dividif(x,y));
disp(dividifHermite(x,y))

function yq = Lagrange(x,y,xq)
    n = length(x);
    if n ~= size(y), error("dati inconsistenti"); end
    %DA FARE
    controllo(xq)
    %-------
    yq = zeros(length(xq));
    for i = 1:n
        yq = yq + y(i)*Lin(x,xq,i);
    end
    return
end

function L = Lin(x,xq,i)
    n = length(x)-1;
    xi = x(i);
    x = x([1:i-1,i+1:n+1]);
    L = ones(length(xq));
    for j=1:n
        L = L*(xq-x(j))/(xi-x(j));
    end
    return
end

function yq = Newton(x,y,xq)
    if length(x) ~= length(y), error("dati inconsistenti"); end
    df = dividif(x,y);
    n = length(df);
    yq = ones(length(xq))*df(n);
    for i = n-1:-1:1
        yq = yq*(xq-x(i))+df(i);
    end
    return
end

function df = dividif(x,y)
    df = y;
    n = length(x);
    for i = 1:n-1
        for j = n:-1:i+1
            df(j) = (df(j)-df(j-1))/(x(j)-x(j-i));
        end
    end
    return
end

function yq = Hermite(x,y,xq)
    if length(x) ~= length(y), error("dati inconsistenti"); end
    df = dividifHermite(x,y);
    n = (length(df)-2)/2;
    yq = ones(length(xq))*df(2*n+2);
    for i = 2*n+1:-1:1
        yq = yq*(xq-x(i))+df(i);
    end
    return
end

function df = dividifHermite(x,y)
    df = y;
    n = (length(df)-2)/2;
    %prima colonna
    for i = 2*n+1:-2:3
        df(i) = (df(i)-df(i-2))/(x(i)-x(i-2));
    end
    %colonne successive
    for i = 2:2*n+1
        for j = 2*n+2:-1:i+1
            df(j) = (df(j)-df(j-1))/(x(j)-x(j-i));
        end
    end
    return
end



