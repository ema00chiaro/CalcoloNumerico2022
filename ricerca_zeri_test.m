addpath("./funcs/Es5");
% f = @functionToPass;
% f1 = @functionToPassDerivata;
% disp("newton");
% [x,i] = newton(f,f1,1,10^-3);
% disp(x);
% disp(i);
% [x,i] = newton(f,f1,1,10^-6);
% disp(x);
% disp(i);
% [x,i] = newton(f,f1,1,10^-9);
% disp(x);
% disp(i);
% [x,i] = newton(f,f1,1,10^-12);
% disp(x);
% disp(i);
% disp("secanti");
% [x,i] = secanti(f,1,0.99,10^-3);
% disp(x);
% disp(i);
% [x,i] = secanti(f,1,0.99,10^-6);
% disp(x);
% disp(i);
% [x,i] = secanti(f,1,0.99,10^-9);
% disp(x);
% disp(i);
% [x,i] = secanti(f,1,0.99,10^-12);
% disp(x);
% disp(i);
% disp("steffensen");
% [x,i] = steffensen(f,1,10^-3);
% disp(x);
% disp(i);
% [x,i] = steffensen(f,1,10^-6);
% disp(x);
% disp(i);
% [x,i] = steffensen(f,1,10^-9);
% disp(x);
% disp(i);
% [x,i] = steffensen(f,1,10^-12);
% disp(x);
% disp(i);

f = @functionToPass2;
f1 = @functionToPass2Derivata;
disp("newton2");
[x,i] = newton(f,f1,1,10^-3);
disp(x);
disp(i);
[x,i] = newton(f,f1,1,10^-6);
disp(x);
disp(i);
[x,i] = newton(f,f1,1,10^-9);
disp(x);
disp(i);
[x,i] = newton(f,f1,1,10^-12);
disp(x);
disp(i);
disp("secanti2");
[x,i] = secanti(f,1,0.99,10^-3);
disp(x);
disp(i);
[x,i] = secanti(f,1,0.99,10^-6);
disp(x);
disp(i);
[x,i] = secanti(f,1,0.99,10^-9);
disp(x);
disp(i);
[x,i] = secanti(f,1,0.99,10^-12);
disp(x);
disp(i);
disp("steffensen2");
[x,i] = steffensen(f,1,10^-3);
disp(x);
disp(i);
[x,i] = steffensen(f,1,10^-6);
disp(x);
disp(i);
[x,i] = steffensen(f,1,10^-9);
disp(x);
disp(i);
[x,i] = steffensen(f,1,10^-12);
disp(x);
disp(i);

function y = functionToPass(x)   
    y = x - cos((pi/2)*x);
end

function y = functionToPassDerivata(x)   
    y = 1 + (pi/2)*sin((pi/2)*x);
end

function y = functionToPass2(x)   
    y = (x - cos((pi/2)*x))^3;
end

function y = functionToPass2Derivata(x)   
    y = (3*(x-cos((pi/2)*x))^2)*(1 + (pi/2)*sin((pi/2)*x));
end
