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
    phi = h(1:n-1)./(h(1:n-1)+h(2:n));
    epsilon = h(2:n)./(h(1:n-1)+h(2:n));
%     diag = ones(1,n-1)*2;
%     x = miaTriLU(diag,phi,epsilon,f);
end

function f = getConstantTerms(x,y)
    n = length(x);
    df2 = (y(2:n)-y(1:n-1))./((x(2:n)-x(1:n-1))); % vettore con f[x0x1]...f[x(n-1)xn]
    df3 = (df2(2:n-1)-df2(1:n-2))./((x(3:n)-x(1:n-2))); %vettore con f[x0x1x2]...f[x(n-2)x(n-1)xn]
    f = df3*6;
end

function x = miaTriLU(a,b,c,f)
    n = length(a);
    for i = 1:n-1
        b(i) = b(i)/a(i);
        a(i+1) = a(i+1)-b(i)*c(i);
    end
%     x = f;
%     %Ly = f
%     for i = 1:n-2
%         x(i+1) = f(i+1)-b(i)*x(i);
%     end
%     %Ux = y
%     x(n) = x(n)/a(n);
%     for i = n-2:-1:1
%         x(i) = (x(i)-c(i)*x(i+1))/a(i);
%     end
    return
end