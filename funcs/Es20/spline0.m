function yq = spline0(x,y,xq)
    n = length(x);
    h = x(2:n)-x(1:n-1);
    %calcolo phi ed epsilon saltando il primo elemento
    %per phi e l'ultimo per epsilon in quanto non necessari
    phi = h(2:end-1)./(h(2:end-1)+h(3:end));
    epsilon = h(2:end-1)./(h(1:end-2)+h(2:end-1));
    diag = ones(1,n-2)*2;
    df = 6*getDiffs(x,y);
    m = miaTriLU(diag,phi,epsilon,df);
    
    yq = zeros(size(xq));

    q = zeros(size(h));
    r = zeros(size(h));

    %calcolo i vettori q ed r
    q(1) = (y(2)-y(1))/h(1) - (h(1)/6)*(m(1));
    r(1) = y(1) - (h(1)^2/6)*m(1);
    for i = 3:n-1
         m1 = m(i-2);
         m2 = m(i-1);
        q(i-1) = (y(i)-y(i-1))/h(i-1) - (h(i-1)/6)*(m2-m1);
        r(i-1) = y(i-1) - (h(i-1)^2/6)*m1;
    end
    q(n-1) = (y(n)-y(n-1))/h(n-1) - (h(n-1)/6)*(-m(n-2));
    r(n-1) = y(n-1) - (h(n-1)^2/6)*m(n-2);

    %per ogni valore di xq, cerco l'intervallo di appartenenza
    %e successivamente calcolo s3(xq)
    for i = 1:length(xq)
        k = 0; %intervallo di appertenza di xq(i)
        for j = 2:n
            if xq(i) >= x(j-1) && xq(i) <= x(j)
                k = j-1;
                break
            end
        end
        if k == 0, error("valore di xq : " + xq(i) + ...
                " non appartiene all'inervallo " +x(1) + "-"+x(end));end
        if k == 1
            m1 = 0;
            m2 = m(1);
        elseif k == length(x)-1
            m1 = m(end);
            m2 = 0;
        elseif k ~= 0
            m1 = m(k-1);
            m2 = m(k);
        end
        if k ~= 0
            yq(i) = ((((xq(i)-x(k))^3)*m2)+(((x(k+1)-xq(i))^3)*m1))/(6*h(k)) + q(k)*(xq(i)-x(k)) + r(k);
        end
    end
    return
end

function diffs = getDiffs(x,y)
    n = length(x);
    df2 = (y(2:n)-y(1:n-1))./(x(2:n)-x(1:n-1));
    diffs = (df2(2:end)-df2(1:end-1))./(x(3:n)-x(1:n-2));
    return
end

function m = miaTriLU(a,b,c,m)
    n = length(a);
    for i = 1:n-1
        b(i) = b(i)/a(i);
        a(i+1) = a(i+1)-b(i)*c(i);
        m(i+1) = m(i+1)-b(i)*m(i);
    end
    m(n) = m(n)/a(n);
    for i = n-1:-1:1
        m(i) = (m(i)-c(i)*m(i+1))/a(i);
    end
    return
end