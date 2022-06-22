a = 0;
b = 1;
integi = zeros(1,9);
erron = ones(1,9)*(1/3)*(exp(3)-1);
for n = 1 : 9
    h = (b-a)/n;
    i = (0:n);
    xi = a+i*h;
    fi = funzione(xi);
    w = weights(n);
    integi(n) = h*sum(fi.*w);
    disp(erron(n));
    disp(integi(n));
    erron(n) = abs(erron(n) - integi(n));
end

n = [1:7,9]';
in = integi([1:7,9])';
err = erron([1:7,9])';
disp(table(n,in,err));


function y = funzione(x)    
    y = exp(3*x);
end

function w = weights(n)
    w = zeros(1,n+1);
    for i = 0:n
        a = poly([0:i-1,i+1:n]);
        a = [a./(n+1:-1:1), 0];
        num = polyval(a,(n));
        d = i - [0:i-1,i+1:n];
        den = prod(d);
        w(i+1) = num/den;
    end
end