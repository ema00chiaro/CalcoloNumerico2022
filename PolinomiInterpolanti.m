x = [1,1,1,2,3,3];
y = [11,10,23,22,34,33];
disp(lagrange([3,6,5,4],[1,2,3,4],[7,2,4]));
function yq = lagrange(x,y,xq)
% yq = Lagrange(x,y,xq)
% 
% Calcola nei punti xq il polinomio di Lagrange interpolante i punti (x,y)
% 
% Input:
%     x - le ascisse di interpolazione (x0,x1,...,xn)
%     y - il valore della funzione nelle ascisse di interpolazione (y0,y1,...,yn)
%     xq - le ascisse in cui si vuole calcolare il polinomio interpolante
% Output:
%     yq - il valore del polinomio interpolante calcolato sulle ascisse xq
    n = length(x);
    if n ~= size(y), error("dati inconsistenti"); end
    if containsDuplicates(x), error("le ascisse non " + ...
            "sono distinte fra loro"); end
    yq = zeros(size(xq));
    for i = 1:n
        yq = yq + y(i)*Lin(x,xq,i);
    end
    return
end

function L = Lin(x,xq,i)
% L = Lin(x,xq,i)
% 
% Calcola nei punti di xq il polinomio di base Lagrange 
% di indice i definito sulle ascisse contenute in x
% 
% Input:
%     x - le ascisse sulle quali definire il polinomio
%     xq - i punti in cui si vuole calore il polinomio
%     i - indice del polinomio di base Lagrange
% Output:
%     L - il polinomio di base Lagrange calcolato nei punti xq
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

function yq = newton(x,y,xq)
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

function yq = hermite(x,y,xq)
    %controlli vari
%     if length(x) ~= length(y), error("dati inconsistenti"); end
%     xapp = x(1:2:length(x));
%     if containsDuplicates(xapp), error("le ascisse non " + ...
%             "sono distinte fra loro" + ...
%             " o non sono scritte nell'ordine x0,x0,x1,x1... ecc"); end
%     xapp2 = x(2:2:length(x));
%     if containsDuplicates(xapp2), error("le ascisse non " + ...
%             "sono distinte fra loro o" + ...
%             " non sono scritte nell'ordine x0,x0,x1,x1... ecc"); end
%     if ~isequal(xapp,xapp2) error("Le ascisse raddoppiate non coincidono " + ...
%             " o non sono scritte nell'ordine x0,x0,x1,x1... ecc"); end
    n = (length(x)-2)/2;
    if length(unique(x)) ~= n+1, error('dati inconsistenti'); end
    for i = 1:2:n-1
        if x(i) ~= x(i+1), error('dati non scritti nella maniera opportuna'); end
    end
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



