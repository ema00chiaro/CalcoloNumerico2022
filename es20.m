x = [2,3,11,6,8];
y = [5,6,7,8,12];
xq = (1:5);
disp(spline0(x,y,xq));

function yq = spline0(x,y,xq)
%VEDERE SE RIUSCIAMO AD USARE SOLAMENTE I 4 VETTORI NECESSARI,
%senza usare h,df2,df3
    n = length(x)-1;
    f = getConstantTerms(x,y);
    h = x(2:n+1)-x(1:n);
    phi = h(2:n-1)./(h(2:n-1)+h(3:n));
    epsilon = h(2:n-1)./(h(1:n-2)+h(2:n-1));
    diag = ones(1,n-1)*2;
    m = miaTriLU(diag,phi,epsilon,f);
    % ora bisogna calcolare s3 sfruttando i valori di m
    % m contiene s3^(II)(xi), calcolare -> s3
    % bisogna capire in quale intervallo si trovano i vari elementi di xq e
    % poi calcolare s3 di quell'intervallo nel punto giusto di xq

end

function f = getConstantTerms(x,y)
    n = length(x);
    df2 = (y(2:n)-y(1:n-1))./((x(2:n)-x(1:n-1))); % vettore con f[x0x1]...f[x(n-1)xn]
    df3 = (df2(2:n-1)-df2(1:n-2))./((x(3:n)-x(1:n-2))); %vettore con f[x0x1x2]...f[x(n-2)x(n-1)xn]
    f = df3*6;
end

function m = miaTriLU(a,b,c,f)
    n = length(a);
    for i = 1:n-1
        b(i) = b(i)/a(i);
        a(i+1) = a(i+1)-b(i)*c(i);
    end
    m = f;
    %Ly = f
    for i = 1:n-1
        m(i+1) = f(i+1)-b(i)*m(i);
    end
    %Ux = y
    m(n) = m(n)/a(n);
    for i = n-1:-1:1
        x(i) = (x(i)-c(i)*m(i+1))/a(i);
    end
    return
end