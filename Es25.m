addpath("./funcs/Es24");
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