% x = [1 2 3 4 5 5.5 7 8 9 9.5 10];
% y = [0 0 0 0.5 0.4 1.2 1.2 0.1 0 0.3 0.6];
% xq = [1.05,2.3,5.6];
% disp(spline0(x,y,xq));

a=0;
b=20*pi;
n = 100;
xq = linspace(a,b,100);
equidistanti = linspace(a,b,n+1);
yeq = functionToPass(equidistanti);


plot(equidistanti,yeq);
hold on;
plot(xq,spline0(equidistanti,yeq,xq));
hold off;

function y = functionToPass(x)
%     y = 1./(2*(2*x.^2 - 2*x + 1));
    y = cos(x);
    return
end

function yq = spline0(x,y,xq)
    n = length(x)-1;
    f = getConstantTerms(x,y);
    h = x(2:n+1)-x(1:n);
    phi = h(2:n-1)./(h(2:n-1)+h(3:n));
    epsilon = h(2:n-1)./(h(1:n-2)+h(2:n-1));
    diag = ones(1,n-1)*2;
    m = miaTriLU(diag,phi,epsilon,f);

    yq = zeros(size(xq));
    for i = 1: length(xq)
        k = 0;
        % fare funzione per questo for
        for j = 1:length(x)-1
            if xq(i) >= x(j) && xq(i) <= x(j+1)
                k = j+1; % k prende l'indice più alto degli estremi dell'intervallo
                break
            end
        end
        if k == 0, error("estremi non trovati");end
        %DA RIVEDERE
        if k == 2 % primo intervallo (x0,x1) che in matlab è (x1,x2)
            mk = m(k-1);
            mk1 = 0;
        elseif k == length(x) % intervallo finale
            mk = 0;
            mk1 = m(k-2);
        else
            mk = m(k-1);
            mk1 = m(k-2);
        end
        qi = (y(k)-y(k-1))/h(k-1) - (h(k-1)/6)*(mk-mk1);
        ri = y(k-1) - ((h(k-1)^2)/6)*mk1;
        yq(i) = (((xq(i)-x(k-1))^3)*mk-((x(k)-xq(i))^3)*mk1)/(6*h(k-1)) + qi*(xq(i)-x(k-1)) + ri;

    end
    return
end

function f = getConstantTerms(x,y)
    n = length(x);
    df2 = (y(2:n)-y(1:n-1))./((x(2:n)-x(1:n-1))); % vettore con f[x0x1]...f[x(n-1)xn]
    df3 = (df2(2:n-1)-df2(1:n-2))./((x(3:n)-x(1:n-2))); %vettore con f[x0x1x2]...f[x(n-2)x(n-1)xn]
    f = df3*6;
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