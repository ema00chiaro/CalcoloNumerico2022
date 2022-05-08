function yq = Lagrange(x,y,xq)
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
