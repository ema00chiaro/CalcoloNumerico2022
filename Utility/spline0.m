
function yq = spline0(x,y,xq)
    n = length(x)-1;
    f = getConstantTerms(x,y);
    h = x(2:n+1)-x(1:n);
    phi = h(2:n-1)./(h(2:n-1)+h(3:n));
%     epsilon = h(2:n-1)./(h(1:n-2)+h(2:n-1));
    epsilon =h(3:n)./(h(2:n-1)+h(3:n));
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
        %DA RIVEDERE======================================================
        %perche in questo momento se una x non è in nessun intervallo
        %lo ingora senza fare niente
        if k == 2 % primo intervallo (x0,x1) che in matlab è (x1,x2)
            mk = m(k-1);
            mk1 = 0;
        elseif k == length(x) % intervallo finale
            mk = 0;
            mk1 = m(k-2);
        elseif k ~= 0
            mk = m(k-1);
            mk1 = m(k-2);
        end
        if k ~= 0
            qi = (y(k)-y(k-1))/h(k-1) - (h(k-1)/6)*(mk-mk1);
            ri = y(k-1) - ((h(k-1)^2)/6)*mk1;
            yq(i) = (((xq(i)-x(k-1))^3)*mk-((x(k)-xq(i))^3)*mk1)/(6*h(k-1)) + qi*(xq(i)-x(k-1)) + ri;
        else
            disp("sono chionzo");
        end
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