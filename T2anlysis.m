function T2anlysis(image,u,u_0,sector,te,si,t2star,a)

subplot(121);
imagesc(image); colormap(gray);hold on;
% set(h1,'color','r');
% set(h2,'color','r');
contour(u,[0 0],'r');
contour(u_0,[0 0],'r');
[C,h] = contour(sector,[0 0]);
L = length(C(1,:));
x1 = C(1,2:L);
y1 = C(2,2:L);
fill(x1,y1,'g');

x = 0:20;
r = -1/t2star;
y =  a*exp(r*x);
subplot(122);
plot(te,si,'*',x,y,'b');
grid on;
axis([0 20 0 180]); %ÉèÖÃ×ø±êÖáµÄ·¶Î§
title('T2* decay curve');
xlabel('TE(Millisecond)');
ylabel('SI');
t2star = round(t2star*100);
t2star = t2star/100;
text(14,170,['T2*=',num2str(t2star)],'color','r');
