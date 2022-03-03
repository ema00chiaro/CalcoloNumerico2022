a = [1,1,1;2,10,5;1,1,1;6,7,4];
b=[3;5;9;6];
c = fattQR(a);
disp(c);

function a = fattQR(a)
    [m,n]=size(a);
    for i = 1:n
        %calcolo alpha con il suo segno
        alpha = norm(a(i:m,i));
        if a(i,i) >= 0, alpha = -alpha; end
        %primo elemento di v
        v1 = a(i,i)-alpha;
        a(i+1:m,i) = a(i+1:m,i)/v1; %memorizzo il vettore v cappuccio
        a(i,i) = alpha;
        %prodotto HA -> HA = A-beta*v*vT*A
        beta = -v1/alpha;
        a(i:m,i+1:n) = a(i:m,i+1:n) - (beta*[1;a(i+1:m,i)])*([1;a(i+1:m,i)]'*a(i:m,i+1:n));
    end
end
% 
% function g = QTb(a,b)
%     g = b;
%     [m,n]=size(a);
%     for i = 1:n
%         v=[1;a(i+1:m,i)];
%         beta = 2/(v'*v);
%         g(i:m) = (eye(m-i+1)-beta*v*v')*g(i:m);
%     end
% end
function g = QTb(a,b)
    [m,n]=size(a);
    g = b;
    %Hn.....H1b
    for i = 1:n
        v=[1;a(i+1:m,i)];
        beta = 2/(v'*v);
        g(i:m) = g(i:m) - (beta * v) * (v'*g(i:m));
        %ricorda che Hi = eye(m-i+1) - beta*vi*vi'
    end
end

function x = minsquaresolverQR(a,g)
    [m,n] = size(a);
    x = g(1:n);
    %Rhat x = g1
    for i = n:-1:1
        x(i) = x(i)/a(i,i);
        x(1:i-1) = x(1:i-1) - a(1:i-1,i)*x(i);
    end
end