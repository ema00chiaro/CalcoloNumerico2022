addpath("./Utility");
a = -2;
b = 3;
n = 40;
xq = linspace(a,b,1001);
equidistanti = linspace(a,b,n+1);
yeq = functionToPass(equidistanti);
cheb = chebyshev(a,b,n);
ycheb = functionToPass(cheb);

%solo con hermite
equiHerm = repelem(equidistanti,2);
chebHerm = repelem(cheb,2);

yeqHerm = functionToPassDifferentiate(equiHerm);
for i = 1:2:length(yeqHerm)
    yeqHerm(i)  = functionToPass(equiHerm(i));
end

ychebHerm = functionToPassDifferentiate(chebHerm);
for i = 1:2:length(ychebHerm)
    ychebHerm(i)  = functionToPass(chebHerm(i));
end
%--------------------

plot(xq,hermite(chebHerm,ychebHerm,xq));


function y = functionToPass(x)
    y = 1./(2*(2*x.^2 - 2*x + 1));
    return
end

function y = functionToPassDifferentiate(x)
    y = (1-2*x)./((2*x.^2 - 2*x + 1).^2);
    return
end

function xc = chebyshev(a,b,n)
    teta = (n:-1:0);
    teta = (teta.*2 + 1).*(pi/(2*(n+1)));
    xc = cos(teta).*((b-a)/2) + ((a+b)/2);
    return
end