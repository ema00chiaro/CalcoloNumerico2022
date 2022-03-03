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