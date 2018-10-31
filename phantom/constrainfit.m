function p = constrainfit(x,y,x0,y0,n)
% c = [1 -2 1 -1];  
% x = linspace(-2,4);  
% y = c(1)*x.^3+c(2)*x.^2+c(3)*x+c(4) + randn(1,100);
% figure,
% plot(x,y,'.b-')

%% fit
% hold on
% x0 = 1;
% y0 = 10;
% n = 3;
C(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
    C(:,j) = x'.*C(:,j+1);
end
d = y';

A = [];
b = [];
Aeq = x0.^(n:-1:0);
beq = y0;

p = lsqlin(C,d,A,b,Aeq,beq);
% yhat = polyval(p,x);
% plot(x0,y0,'gx','linewidth',4) 
% % Plot fitted data
% plot(x,yhat,'r','linewidth',2) 
% hold off