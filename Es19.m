disp(Chebyshev(1,2,2))
function xc = Chebyshev(a,b,n)
    if b <= a, error("a deve essere minore di b");end
    if n <= 0, error("il numero ascisse richieste" + ...
            " deve essere positivo");end

    teta = (n:-1:0);
    teta = (teta.*2 + 1).*(pi/(2*(n+1)));
    xc = cos(teta).*((b-a)/2) + ((a+b)/2);
    return
end
