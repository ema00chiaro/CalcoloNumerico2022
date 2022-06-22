function compareSplines()
    xq = linspace(0,1,10001);
    n = 5;
    ns = (5:5:50);
    errs1 = zeros(1,10);
    errs2 = zeros(1,10);
    for i = 1:10 
        x = linspace(0,1,n+1);
        y = f1(x);
        errs1(i) = norm(abs(spline(x,y,xq)-spline0(x,y,xq)),"inf");
        y = f2(x);
        errs2(i) = norm(abs(spline(x,y,xq)-spline0(x,y,xq)),"inf");
        n = n+5;
    end
    
    disp(table(ns',errs1',errs2'));
end

function y = f1(x)
    y = sin(2*pi*x);
end

function y = f2(x)
    y = cos(2*pi*x);
end