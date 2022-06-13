tols = [10^-2;10^-3;10^-4;10^-5;10^-6];
ris = zeros(1,5);
vals = zeros(1,5);
errs = zeros(1,5);
for i = 1:9
    if i == 8
    else
        disp("ITERAZIONE-------- " + i);
        [ris(1),errs(1),vals(1)] = composita(@functionToPass,0,1,i,10^-2);
        [ris(2),errs(2),vals(2)] = composita(@functionToPass,0,1,i,10^-3);
        [ris(3),errs(3),vals(3)] = composita(@functionToPass,0,1,i,10^-4);
        [ris(4),errs(4),vals(4)] = composita(@functionToPass,0,1,i,10^-5);
        [ris(5),errs(5),vals(5)] = composita(@functionToPass,0,1,i,10^-6);
        disp(table(tols,ris',errs',vals'));
    end
end




function y =  functionToPass(x)
    y = sin(1./(0.1+x));
end

function [If,err,nfeval] = composita(fun,a,b,n,tol)
    err = tol + 1;
    nfeval = 0;
    m = n;
    If = 0;
    while (err > tol)
        m = 2*m;
        h = (b-a)/m;
        i = (0:m);
        xi = a+i*h;
        if m == 2*n % first time
            fi = fun(xi);
            nfeval = nfeval + m+1;
        else
            fi = repelem(fi,2);
            fi = fi(1:m+1);
            fi(2:2:m) = fun(xi(2:2:m));
            nfeval = nfeval + m/2;
        end
        if mod(n,2) == 0
            mu = 2;
        else
            mu = 1;
        end

        w = weights(n);
        Ieven = 0;

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
%         disp("-----------");
%         disp("m = " + m);
%         disp("If = " + If);
%         disp("errore " + err);
%         disp("nfeval " + nfeval);
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