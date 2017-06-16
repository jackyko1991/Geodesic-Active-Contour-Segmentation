function T2curve(t2star,te,si,a)
x = 0:20;
r = -1/t2star;
y =  a*exp(r*x);
plot(te,si,'*',x,y,'b');
grid on;
axis([0 20 0 500]); %ÉèÖÃ×ø±êÖáµÄ·¶Î§
title('T2* decay curve','FontSize',14);
xlabel('TE(Millisecond)','FontSize',14,'color','b');
ylabel('SI','FontSize',14,'color','b');
t2star = round(t2star*100);
t2star = t2star/100;
text(14,180,['T2*=',num2str(t2star)],'FontSize',14,'color','r');
