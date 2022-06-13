composita(@functionToPass,0,pi,2,10^-5);

function y =  functionToPass(x)
    y = sin(1./(0.1+x));
end

function [If,err,nfeval] = composita(fun,a,b,n,tol)
    En = tol + 1;
    nfeval = 0;
    m = n;
    ite = 0;
    while (En > tol)
        m = 2*m;
        h = (b-a)/m;
        i = (0:m);
        xi = a+i*h;
        fi = fun(xi);
        nfeval = nfeval + m+1;

        if mod(n,2) == 0
            mu = 2;
        else
            mu = 1;
        end
        %ora ci calcoliamo Inm che sarebbe Ikn
        w = weights(n);
        Inm = 0;
        Inm2 = 0;
%         index = 1;
%         for k = 1:1:(m/2)/n
%             adesso = (xi((n*k)*2 +1)-xi((n*k-n)*2 +1))/n;
%             somma = 0;
%             for j = 1:1:n+1
%                 somma = somma + (w(j)*fi(index));
%                 index = index + 2;
%             end
%             adesso = adesso * somma;
%             Inm2 = Inm2 + adesso;
%             index = index - 2;
%         end
% 
%         index = 1;
%         for k = 1:1:m/n
%             adesso = ((xi(n*k +1)-xi(n*k-n +1))/n);
%             somma = 0;
%             for j = 1:1:n+1
%                 somma = somma + (w(j)*fi(index));
%                 index = index + 1;
%             end
%             adesso = adesso * somma;
%             Inm = Inm + adesso;
%             index = index - 1;
%         end
    
        k = 1;
        for j = 1:2:m+1
%             disp("peso: " + w(k));
%             disp("fi: " + fi(j));

            Inm2 = Inm2 + fi(j)*w(k);
            k = k+1;
            if k-1 > n
                currentK = k;
                k = 1;
                if j ~=  m+1
                    Inm2 = Inm2 + fi(j)*w(currentK-1);
                    k = 2;
                end
            end
        end
        Inm2 = Inm2 *(b-a)/(m/2);
        
        k = 1;
        for j = 1:1:m+1
%             disp("peso: " + w(k));
%             disp("fi: " + fi(j));
            Inm = Inm + fi(j)*w(k);
            k = k+1;
            if k-1 > n
                currentK = k;
                k = 1;
                if j ~=  m+1
                    Inm = Inm + fi(j)*w(currentK-1);
                    k = 2;
                end
            end
        end
        Inm = Inm * ((b-a)/(m));


        En = abs(Inm-Inm2)/(2^(n+mu)-1);
        ite = ite +1;
        disp("iteraz-----------" + ite);
        disp("m = " + m);
        disp("Inm = " + Inm);
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