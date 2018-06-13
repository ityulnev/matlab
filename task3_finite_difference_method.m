clear all
%Finite Difference Method

%Meshing
step=0.1;
T=10;
L=20;

t=0:step:T;
x=0:step:L;
dt=t(2)-t(1);
dx=x(2)-x(1);

c=1;%3e8;  %Light Speed
C=c*dt/dx; %Courant Number

%Init Matrixes
[xx,tt]=meshgrid(x,t);
u_exact=zeros(length(t),length(x));


%Tasks from http://hplgit.github.io/INF5620/doc/pub/sphinx-wave/._main_wave001.html
task='A';
switch task
    case 'A'
    %A) u_exact=x(L-x)(1+t/2) 
    I=x.*(L-x);%Initial Signal at t=0
    V=0.5.*x.*(L-x);%d/dt(u_e)
    f=(c^2)*2.*(1+0.5.*tt);   
    u_exact=xx.*(L-xx).*(1+tt./2);
     
    case 'B'%Correct 
    %B) u_exact=sin(pi*x/L)*cos(pi*c*t/L)
    I=sin(pi.*x./L);%Initial Signal at t=0
    V=x.*0;
    f=tt.*0;    
    u_exact=sin(pi.*xx./L).*cos(pi.*c.*tt./L);
  
    case 'C'%result is exactly the negative of the u_exact
    %C) u_exact=x(L-x)sin(c*pi*t/L)
    I=x.*0;%Initial Signal at t=0
    V=x.*(L-x).*pi*c/L;
    f=((c^2)*2-xx.*(L-xx).*((pi*c/L)^2)).*sin(c*pi.*tt./L);    
    u_exact=xx.*(L-xx).*sin(c*pi.*tt./L);
     
    case 'C2'
    %C) u_exact=x(L-x)sin(t)
    I=x.*0;%Initial Signal at t=0
    V=x.*(L-x);
    f=(c^2)*2-xx.*(L-xx).*sin(tt);    
    u_exact=xx.*(L-xx).*sin(tt);    
end


%u=doSplitStep(I,x,t,C);
u=doGeneralSplitStep(I,V,f,x,t,C);
%u_exact=doSplitStep(I,x,t,C);    
 
%Plot
figure;
plot1=surf(x,t,u);
hold on
plot2=surf(x,t,u_exact);
hold off
xlabel('x')
ylabel('t')
set(plot1,'LineStyle','none')
set(plot2,'LineStyle','none','facealpha',0.7,'facecolor','black')
legend('Finite Diff Method','Exact Solution')


