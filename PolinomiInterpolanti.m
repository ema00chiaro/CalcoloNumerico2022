x = [1,1,1,2,3,3];
y = [11,10,23,22,34,33];
disp(Lagrange([3,4,5],[1,2,3],[7,2,4]));
function yq = Lagrange(x,y,xq)
    n = length(x);
    if n ~= size(y), error("dati inconsistenti"); end
    if containsDuplicates(xq), error("le ascisse non " + ...
            "sono distinte fra loro"); end
    yq = zeros(size(xq));
    for i = 1:n
        yq = yq + y(i)*Lin(x,xq,i);
    end
    return
end

function L = Lin(x,xq,i)
    n = length(x)-1;
    xi = x(i);
    x = x([1:i-1,i+1:n+1]);
    L = ones(size(xq));
    for j=1:n
        L = L.*(xq-x(j))/(xi-x(j));
    end
    return
end

function containsDuplicates = containsDuplicates(x)
    containsDuplicates = length(x) ~= length(unique(x));
    return
end

function yq = Newton(x,y,xq)
    if length(x) ~= length(y), error("dati inconsistenti"); end
    if containsDuplicates(xq), error("le ascisse non " + ...
            "sono distinte fra loro"); end
    df = dividif(x,y);
    n = length(df);
    yq = ones(size(xq))*df(n);
    for i = n-1:-1:1
        yq = yq.*(xq-x(i))+df(i);
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
    %controlli vari
    if length(x) ~= length(y), error("dati inconsistenti"); end
    xapp = x(1:2:length(x));
    if containsDuplicates(xapp), error("le ascisse non " + ...
            "sono distinte fra loro" + ...
            " o non sono scritte nell'ordine x0,x0,x1,x1... ecc"); end
    xapp2 = x(2:2:length(x));
    if containsDuplicates(xapp2), error("le ascisse non " + ...
            "sono distinte fra loro o" + ...
            " non sono scritte nell'ordine x0,x0,x1,x1... ecc"); end
    if ~isequal(xapp,xapp2) error("Le ascisse raddoppiate non coincidono " + ...
            " o non sono scritte nell'ordine x0,x0,x1,x1... ecc"); end
    %fine controlli
    df = dividifHermite(x,y);
    n = (length(df)-2)/2;
    yq = ones(size(xq))*df(2*n+2);
    for i = 2*n+1:-1:1
        yq = yq.*(xq-x(i))+df(i);
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



