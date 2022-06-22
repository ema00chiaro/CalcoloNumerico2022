function yq = hermite(x,y,y1,xq)
% yq = hermite(x,y,y1,xq)
%
% Funzione che calcola il polinomio interpolante di Hermite nei
% punti xq.
% Input:
% 	x - vettore che contiene le ascisse di interpolazione
%	y - valori che assume la funzione nei punti x
%	y1 - valori che assume la derivata della funzione nei punti x
%	xq - valori dove vogliamo calcolare il polinomio
% 		interpolante
% Output:
%	yq - valori che assume il polinomio interpolante in
% 		 corrispondenza di xq.

    n = length(x)-1;
    if length(x) ~= length(y) || length(x) ~= length(y1)
         error('dati inconsistenti');
    end
    if length(unique(x)) ~= n+1, error('dati inconsistenti'); end
    
    x = repelem(x,2);
    y = repelem(y,2);
    y(2:2:end) = y1(1:end);
    df = dividifHermite(x,y);
    yq = ones(size(xq))*df(2*n+2);
    for i = 2*n+1:-1:1
        yq = yq.*(xq-x(i))+df(i);
    end
    return
end

function df = dividifHermite(x,y)
%df = dividifHermite(x,y)
%
% Funzione che calcola il vettore delle differenze divise per il %polinomio interpolante di Hermite.
% Input:
%	x - valori delle ascisse di interpolazione
%	y - valori della funzione in corrispondenza di x
% Output:
%	df - vettore contenente le differenze divise

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
