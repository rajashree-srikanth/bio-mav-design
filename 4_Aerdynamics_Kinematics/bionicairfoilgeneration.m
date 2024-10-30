S = [3.9362 -0.7705 0.8485 1];
A = [-29.4861 66.4565 -59.8060 19.0439];
zc = zeros(5);
zt = zeros(5);
zct = 0;
ztt = 0;
zcp = 0;
ztp = 0;
zu = 0;
zl = 0;
xcord = [];
yl = [];
yu = [];

for y = 0:0.0414:0.414
    zcmax = (0.18)/(1 + (7.31)*(y^2.77));
    ztmax = (0.1)/(1 + (14.86)*(y^3.52));
    x = 0;
    for x = 0:0.01:1
        xcord(end+1) = x;
            for n = 1:4
            zcp = zcmax*(x)*(1-x)*(S(n)*((2*x - 1)^(n-1)));
            zc(n) = zcp;
            ztp = ztmax*(A(n)*(x^(n + 1) - x^0.5));
            zt(n) = ztp;
            if n<=3
                zct = zct + zc(n);
            end    
            ztt = ztt + zt(n);
        end
%         disp(zct);
%         disp(ztt);
        zu = zct + ztt;
        zl = zct - ztt;
%         disp("Upper coordinate at x = "); 
%         disp(x);
%         disp(zu);
%         disp("Lower coordinate at x = ");
%         disp(x);
%         disp(zl);
        yl(end+1) = zl;
        yu(end+1) = zu;
       
    end
figure(1)
plot(xcord, yl)
hold on
%figure(2)
plot(xcord, yu)
end
        
        