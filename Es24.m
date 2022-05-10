% n = 2; % n = [(1:7),9];

for n = 1:7
    disp("------------- "+ n +" ------------------");
    disp(weights(n));
end
disp("------------- "+ 9 + " ------------------");
disp(weights(9));

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