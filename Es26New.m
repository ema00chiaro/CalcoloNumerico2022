composita(@functionToPass,0,1,7,10^-6);

function y =  functionToPass(x)
    y = exp(3.*x);
end

function [If,err,nfeval] = composita(fun,a,b,n,tol)
    err = tol + 1;
    nfeval = 0;
    m = n;
    ite = 0;
    If = 0;
    while (err > tol)
        m = 2*m;
        h = (b-a)/m;
        i = (0:m);
        xi = a+i*h;
        if m == 2*n,
            fi = fun(xi);
        else
            fi = 
        end
        nfeval = nfeval + m+1;

        if mod(n,2) == 0
            mu = 2;
        else
            mu = 1;
        end
        %ora ci calcoliamo Inm che sarebbe Ikn
        w = weights(n);
        Ieven = 0; %lui bisogna usare quelli pari + lo zero

        i = 1;
        k = 1;
        while i <= m+1
            Ieven = Ieven + w(k)*fi(i);
            i = i+2;
            k = k+1;
            if k == n+2
                k = 1;
                i = i-2;
            end
        end
        Ieven = Ieven * (b-a)/(m/2);
        
        i = 1;
        k = 1;
        while i <= m+1
            If = If + w(k)*fi(i);
            i = i+1;
            k = k+1;
            if k == n+2
                k = 1;
                i = i-1;
            end
        end
        If = If * (b-a)/(m);
    
        err = abs(Ieven-If)/(2^(n+mu)-1);
        ite = ite +1;
        disp("iteraz-----------" + ite);
        disp("m = " + m);
        disp("If = " + If);
        disp("errore " + En);
    end
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