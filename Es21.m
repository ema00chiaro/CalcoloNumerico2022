addpath("./Utility");
a = -2;
b = 3;
n = 32;
xq = linspace(a,b,1001);
equidistanti = linspace(a,b,n+1);
yeq = functionToPass(equidistanti);
cheb = chebyshev(a,b,n);
ycheb = functionToPass(cheb);
v = [4,8,16,32,40];

% yqe = lagrange(equidistanti,yeq,xq);
% yqc = lagrange(cheb,ycheb,xq);

% for j = 1:5
%     n = v(j);
%     disp("----------------" + n)
%     equidistanti = linspace(a,b,n+1);
%     yeq = functionToPass(equidistanti);
%     cheb = chebyshev(a,b,n);
%     ycheb = functionToPass(cheb);
% 
%     disp("lagrange");
%     disp(interpolErrLagrange(equidistanti,yeq,xq,@functionToPass));
%     disp(interpolErrLagrange(cheb,ycheb,xq,@functionToPass));
%     disp("newton");
%     disp(interpolErrNewton(equidistanti,yeq,xq,@functionToPass));
%     disp(interpolErrNewton(cheb,ycheb,xq,@functionToPass));
%     
%     
%     equiHerm = repelem(equidistanti,2);
%     chebHerm = repelem(cheb,2);
%     
%     yeqHerm = functionToPassDifferentiate(equiHerm);
%     for i = 1:2:length(yeqHerm)
%         yeqHerm(i)  = functionToPass(equiHerm(i));
%     end
%     
%     ychebHerm = functionToPassDifferentiate(chebHerm);
%     for i = 1:2:length(ychebHerm)
%         ychebHerm(i)  = functionToPass(chebHerm(i));
%     end
%     
%     disp("hermite");
%     disp(interpolErrHermite(equiHerm,yeqHerm,xq,@functionToPass));
%     disp(interpolErrHermite(chebHerm,ychebHerm,xq,@functionToPass));
%     
%     
%     %DA RIVEDERE
%     disp("slpine0");
%     disp(interpolErrSpline0(equidistanti,yeq,xq,@functionToPass));
%     disp(interpolErrSpline0(cheb,ycheb,xq,@functionToPass));
%     %-----------
%     
%     disp("spline");
%     disp(interpolErrSpline(equidistanti,yeq,xq,@functionToPass));
%     disp(interpolErrSpline(cheb,ycheb,xq,@functionToPass));
% end

eel = zeros(1,5);
ecl = zeros(1,5);
een = zeros(1,5);
ecn = zeros(1,5);
eeh = zeros(1,5);
ech = zeros(1,5);
ees0 = zeros(1,5);
ecs0 = zeros(1,5);
ees = zeros(1,5);
ecs = zeros(1,5);

for j = 1:5
    n = v(j);
    disp("----------------" + n)
    equidistanti = linspace(a,b,n+1);
    yeq = functionToPass(equidistanti);
    cheb = chebyshev(a,b,n);
    ycheb = functionToPass(cheb);

    disp("lagrange");
    eel(j)=interpolErrLagrange(equidistanti,yeq,xq,@functionToPass);
    ecl(j)=interpolErrLagrange(cheb,ycheb,xq,@functionToPass);
    disp("newton");
    een(j)=interpolErrNewton(equidistanti,yeq,xq,@functionToPass);
    ecn(j)=interpolErrNewton(cheb,ycheb,xq,@functionToPass);
    
    
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
    
    disp("hermite");
    eeh(j)=interpolErrHermite(equiHerm,yeqHerm,xq,@functionToPass);
    ech(j)=interpolErrHermite(chebHerm,ychebHerm,xq,@functionToPass);
    
    
    %DA RIVEDERE
    disp("slpine0");
    ees0(j)=interpolErrSpline0(equidistanti,yeq,xq,@functionToPass);
    ecs0(j)=interpolErrSpline0(cheb,ycheb,xq,@functionToPass);
    %-----------
    
    disp("spline");
    ees(j)=interpolErrSpline(equidistanti,yeq,xq,@functionToPass);
    ecs(j)=interpolErrSpline(cheb,ycheb,xq,@functionToPass);
end

n = v';
eel = eel';
ecl = ecl';
een = een';
ecn = ecn';
eeh = eeh';
ech = ech';
ees0 = ees0';
ecs0 = ecs0';
ees = ees';
ecs = ecs';
format long e;
disp(table(n,eel,een,eeh,ees0,ees));
disp(table(n,ecl,ecn,ech,ecs0,ecs));

% plot(equidistanti,yeq);
% hold on;
% plot(xq,spline0(equidistanti,yeq,xq));
% hold off;

function y = functionToPass(x)
    y = 1./(2*(2*x.^2 - 2*x + 1));
    return
end

function y = functionToPassDifferentiate(x)
    y = (1-2*x)/(2*x.^2 - 2*x + 1).^2;
    return
end

function xc = chebyshev(a,b,n)
    teta = (n:-1:0);
    teta = (teta.*2 + 1).*(pi/(2*(n+1)));
    xc = cos(teta).*((b-a)/2) + ((a+b)/2);
    return
end

function e = interpolErrLagrange(x,y,xq,f)
    e = 1 + lebesgue(x,xq);
    e = e * norm(abs(feval(f,xq)-lagrange(x,y,xq)),"inf");
%     e = norm(abs(feval(f,xq)-lagrange(x,y,xq)),"inf");
end

function e = interpolErrNewton(x,y,xq,f)
    e = 1 + lebesgue(x,xq);
    e = e * norm(abs(feval(f,xq)-newton(x,y,xq)),"inf");
%     e = norm(abs(feval(f,xq)-newton(x,y,xq)),"inf");
end

function e = interpolErrHermite(x,y,xq,f)
    e = 1 + lebesgue(unique(x),xq);
    e = e * norm(abs(feval(f,xq)-hermite(x,y,xq)),"inf");
%     e = norm(abs(feval(f,xq)-hermite(x,y,xq)),"inf");
end

function e = interpolErrSpline0(x,y,xq,f)
    e = 1 + lebesgue(x,xq);
    e = e * norm(abs(feval(f,xq)-spline0(x,y,xq)),"inf");
%     e = norm(abs(feval(f,xq)-spline0(x,y,xq)),"inf");
end

function e = interpolErrSpline(x,y,xq,f)
    e = 1 + lebesgue(x,xq);
    e = e * norm(abs(feval(f,xq)-spline(x,y,xq)),"inf");
%     e = norm(abs(feval(f,xq)-spline(x,y,xq)),"inf");
end

function cLebesgue = lebesgue(x,xq)
    cLebesgue = zeros(size(xq));
    for i = 1:length(x)
        cLebesgue = cLebesgue + abs(Lin(x,xq,i));
    end
    cLebesgue = norm(cLebesgue,"inf");
    return
end

function L = Lin(x,xq,i)
% L = Lin(x,xq,i)
% 
% Calcola nei punti di xq il polinomio di base Lagrange 
% di indice i definito sulle ascisse contenute in x
% 
% Input:
%     x - le ascisse sulle quali definire il polinomio
%     xq - i punti in cui si vuole calore il polinomio
%     i - indice del polinomio di base Lagrange
% Output:
%     L - il polinomio di base Lagrange calcolato nei punti xq
    n = length(x)-1;
    xi = x(i);
    x = x([1:i-1,i+1:n+1]);
    L = ones(size(xq));
    for j=1:n
        L = L.*(xq-x(j))/(xi-x(j));
    end
    return
end