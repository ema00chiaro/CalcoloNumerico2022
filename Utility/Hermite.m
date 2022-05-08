function yq = Hermite(x,y,xq)
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
    if length(unique(x)) ~= n, error('dati inconsistenti'); end
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