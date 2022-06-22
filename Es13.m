addpath("./funcs/Es12");
addpath("./funcs","linsis.m");
% a = [1,1,1;2,10,5;1,1,1;6,7,4];
% b=[3;5;9;6];
% [x,nr] = miaqr(a,b);
% disp(x);
% disp(nr);

A = [ 1 3 2; 3 5 4; 5 7 6; 3 6 4; 1 4 2 ];
b = [ 15 28 41 33 22 ]';
% D = diag(1:5);
% D = eye(5,5)+[2,0,0,0,0;0,0,0,0,0;0,0,0,0,0;0,0,0,0,0;0,0,0,0,2];
D = -eye(5,5)*6;
% disp(A\b);
[x,nr] = miaqr(A,b);
disp(x);
disp("err: " + abs(A\b-x)./abs(A\b));
disp(nr);
% disp(D*A\D*b);
[x,nr] = miaqr(D*A,D*b);
disp(x);
disp("err: " + abs((D*A)\(D*b)-x)./abs((D*A)\(D*b)));
disp(nr);


% A = rand(7,5)*10;
% b = rand(7,1)*10;
% sol = miaqr(A,b);
% writematrix(A,'matrici.csv','WriteMode','append')
% writematrix("-",'matrici.csv','WriteMode','append')
% writematrix(b,'matrici.csv','WriteMode','append')
% writematrix("-",'matrici.csv','WriteMode','append')
% writematrix(sol,'matrici.csv','WriteMode','append')
% writematrix("-",'matrici.csv','WriteMode','append')
% % A = rand(7,5)*10;
% % b = rand(7,1);
% disp(A\b);
% disp(miaqr(A,b));


