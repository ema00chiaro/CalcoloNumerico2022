disp(Chebyshev(2,1,2))
function xc = Chebyshev(n,a,b)
    teta = (n:-1:0);
    teta = (teta.*2 + 1).*(pi/(2*(n+1)));
    xc = cos(teta).*((b-a)/2) + ((a+b)/2);
    return
end
