clear all;

A=[1 -1 1
1 0 0
1 1 1];

b=[1/26 1 1/26]';

a=A\b

x=[-1.0:0.01:1.0];

exact=1./(1+25*x.^2);

u=a(1) +a(2)*x + a(3)*(x.^2);


A=[1 -1 1 -1 1 -1 1
1 -2/3 (-2/3)^2 (-2/3)^3  (-2/3)^4  (-2/3)^5  (-2/3)^6 
1 -1/3 (-1/3)^2 (-1/3)^3  (-1/3)^4  (-1/3)^5  (-1/3)^6 
1 0 0 0 0 0 0 
1 1/3 (1/3)^2 (1/3)^3  (1/3)^4  (1/3)^5  (1/3)^6
1 2/3 (2/3)^2 (2/3)^3  (2/3)^4  (2/3)^5  (2/3)^6 
1 1 1 1 1 1 1
];

b=[1/26 1/(1+25*(-2/3)^2) 1/(1+25*(-1/3)^2) 1 1/(1+25*(1/3)^2) 1/(1+25*(2/3)^2) 1/26]';


a=A\b;

u2=a(1) +a(2)*x + a(3)*(x.^2) + a(4)*(x.^3) + a(5)*(x.^4) + a(6)*(x.^5) + a(7)*(x.^6);

A=[1 -1 1 -1 1 -1 1 -1 1 -1 1
1 -.2 (-.2)^2 (-.2)^3 (-.2)^4 (-.2)^5 (-.2)^6 (-.2)^7 (-.2)^8 (-.2)^9 (-.2)^10 
1 -.4 (-.4)^2 (-.4)^3 (-.4)^4 (-.4)^5 (-.4)^6 (-.4)^7 (-.4)^8 (-.4)^9 (-.4)^10
1 -.6 (-.6)^2 (-.6)^3 (-.6)^4 (-.6)^5 (-.6)^6 (-.6)^7 (-.6)^8 (-.6)^9 (-.6)^10
1 -.8 (-.8)^2 (-.8)^3 (-.8)^4 (-.8)^5 (-.8)^6 (-.8)^7 (-.8)^8 (-.8)^9 (-.8)^10 
1 0 0 0 0 0 0 0 0 0 0 
1 .8 (.8)^2 (.8)^3 (.8)^4 (.8)^5 (.8)^6 (.8)^7 (.8)^8 (.8)^9 (.8)^10 
1 .6 (.6)^2 (.6)^3 (.6)^4 (.6)^5 (.6)^6 (.6)^7 (.6)^8 (.6)^9 (.6)^10
1 .4 (.4)^2 (.4)^3 (.4)^4 (.4)^5 (.4)^6 (.4)^7 (.4)^8 (.4)^9 (.4)^10
1 .2 (.2)^2 (.2)^3 (.2)^4 (.2)^5 (.2)^6 (.2)^7 (.2)^8 (.2)^9 (.2)^10 
1 1 1 1 1 1 1 1 1 1 1
];

b=[1/26  1/(1+25*(-.2)^2) 1/(1+25*(-.4)^2) 1/(1+25*(-.6)^2) ...
1/(1+25*(-.8)^2) 1 1/(1+25*(.8)^2) 1/(1+25*(.6)^2) ...
1/(1+25*(.4)^2) 1/(1+25*(.2)^2) 1/26]';

a=A\b;

u3=a(1) +a(2)*x + a(3)*(x.^2) + a(4)*(x.^3) ...
+ a(5)*(x.^4) + a(6)*(x.^5) + a(7)*(x.^6) ...
+ a(8)*(x.^7) + a(9)*(x.^8) + a(10)*(x.^9) ...
+ a(11)*(x.^10) ;




figure (1)

plot(x,exact, '-k', x,u,'--k',x,u2,':k',x,u3,'-.k','LineWidth',2);
axis([-1 1 -0.5 2]);
xlabel('x');
ylabel('u(x)')
legend('Actual function','2nd order polynomial','6th order polynomial'...
  ,'10th order polynomial','Location','North');
set(gca,'fontsize',18);
%set(gca,'units','points');
%set(gca,'position',[58 60 400 300]);

print -depsc -r150 runge.eps

%figure(2)
%plot(x,up,'-');
