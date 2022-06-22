addpath("./funcs/Es8");
addpath("./funcs","linsis.m");
% l = [1,0,0,0;1/2,1,0,0;3/4,1/2,1,0;0,1/2,1/2,1];
% u=[4,0,1,1;0,2,7/2,1/2;0,0,1/2,0;0,0,0,-1/4];
% a=[4,0,1,1;3,1,3,1;0,1,2,0;2,2,4,1];
% p=[1,0,0,0;0,0,0,1;0,1,0,0;0,0,1,0];
% b=[3;5;9;6];
% [a,pos] = fattPALU(a);
% disp(a);
% x = solvePALU(a,b,pos);
% disp(x);

A = rand(5,5)*100;
b = rand(5,1)*10;
sol = mialu(A,b);
writematrix(A,'matrici.csv')
writematrix("-",'matrici.csv','WriteMode','append')
writematrix(b,'matrici.csv','WriteMode','append')
writematrix("-",'matrici.csv','WriteMode','append')
writematrix(sol,'matrici.csv','WriteMode','append')
writematrix("-",'matrici.csv','WriteMode','append')
% disp(abs(A\b - mialu(A,b))./abs(A\b));

% [A,b] = linsis(10,1);
% disp(cond(A));
% sol = mialu(A,b);
% disp(sol);
% writematrix(A,'matrici.csv')
% writematrix("-",'matrici.csv','WriteMode','append')
% writematrix(b,'matrici.csv','WriteMode','append')
% writematrix("-",'matrici.csv','WriteMode','append')
% writematrix(sol,'matrici.csv','WriteMode','append')
% writematrix("-",'matrici.csv','WriteMode','append')
% [A,b] = linsis(10,10);
% disp(cond(A));
% sol = mialu(A,b);
% disp(sol);
% x = sol;
% disp(table(A,b,x))
% writematrix(A,'matrici.csv','WriteMode','append')
% writematrix("-",'matrici.csv','WriteMode','append')
% writematrix(b,'matrici.csv','WriteMode','append')
% writematrix("-",'matrici.csv','WriteMode','append')
% writematrix(sol,'matrici.csv','WriteMode','append')
% writematrix("-",'matrici.csv','WriteMode','append')