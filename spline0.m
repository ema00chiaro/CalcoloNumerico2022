function yq = spline0(x,y,xq)
    n = length(x);
    h = x(2:n-1) - x(1:n-2); % h = xi - x(i-1) con i = 1 ... n-1 e x = x0 ... xn
    df = (y(2:n)-y(1:n-1))./((x(2:n)-x(1:n-1))); % vettore con f[x0x1]...f[x(n-1)xn]
    terminenoto = (df(2:n-1)-df(1:n-2)).*6; %termine noto 6*(f[x(i-1)xi]-f[xi,x(i+1)])
end