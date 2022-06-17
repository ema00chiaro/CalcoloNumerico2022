function yq = spline02(x,y,xq)
    n = length(x);
    h = x(2:n)-x(1:n-1);
    phi = h(1:end-1)./(h(1:end-1)+h(2:end));
    phi = phi(2:end);
    epsilon = h(2:end)./(h(1:end-1)+h(2:end));
    epsilon = epsilon(1:end-1);
    diag = ones(1,n-2)*2;
    f = 6*getDiffs(x,y);
    m = miaTriLU(diag,phi,epsilon,f);
    
    yq = zeros(size(xq));
    for i = 1:length(xq)
        k = 0; %intervallo di appertenza di xq(i)
        for j = 2:n
            if xq(i) >= x(j-1) && xq(i) <= x(j)
                k = j-1;
                break
            end
        end
        m1 = 0;
        m2 = 0;
        if k == 1
            m1 = 0;
            m2 = m(1);
        elseif k == length(x)-1
            m1 = m(end);
            m2 = 0;
        else
            m1 = m(k-1);
            m2 = m(k);
        end
        qj = (y(k+1)-y(k))/h(k) - (h(k)/6)*(m2*m1);
        rj = y(k) - (h(k)^2/6)*m1;
        yq(i) = ((((xq(i)-x(k))^3)*m2)+(((x(k+1)-xq(i))^3)*m1))/(6*h(k)) + qj*(xq(i)-x(k)) + rj;
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