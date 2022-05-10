n = 1e+4;
xi = (0:n)./n;
yi = functionToPassPerturb(xi);
a = 0;
b = 1;
% m = 15;
% p = polyfit(xi,yi,m);
% xq = linspace(a,b,1000);
e = zeros(1,15);
for i = 1:15
    p = polyfit(xi,yi,i);
    e(i) = norm(functionToPass(xi)-polyval(p,xi),"inf");
end
m = (1:15);
semilogy(m,e);
xlabel("m");
ylabel("e = || f-p ||");


% plot(xq,polyval(p,xq));
% hold on;
% plot(xq,functionToPass(xq));
% hold off;

function y = functionToPass(x)
    y = sin(pi*x.^2);
end

function y = functionToPassPerturb(x)
y = functionToPass(x) + 1e-1 * rand(size(x));
end