function yq = Newton(x,y,xq)
    if length(x) ~= length(y), error("dati inconsistenti"); end
    if containsDuplicates(x), error("le ascisse non " + ...
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

function containsDuplicates = containsDuplicates(x)
% containsDuplicates = containsDuplicates(x)
% 
% Funzione che verifica se un vettore contiene elementi duplicati
% Input:
%     x - il vettore di partenza
% Output:
%     containsDuplicates - se assume valore 1 rappresenta il fatto che il vettore contenga 
%                          o meno elementi duplicati, altrimenti assumerà valore 0
    containsDuplicates = length(x) ~= length(unique(x));
    return
end