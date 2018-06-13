function fnct=doGeneralSplitStep(I0,V0,f,x,t,C)
%Does the Split Step method for DFE u_tt=c^2*u_xx+f(t,x)
%I0 is x_i at t=0, f(t,x) is the Source term and V(x) is d/dt(u) with t->0
X=length(x);
T=length(t);
dt=t(2)-t(1);

hist=zeros(T,X);
u=zeros(1,X);
unext=zeros(1,X);

%Initial Conditions at t=0 and t=1 (corresponds to t(1) and t(2))
uprev=I0;
u(2:X-1)=uprev(2:X-1)-dt.*V0(2:X-1)+0.5*(C^2).*(uprev(3:X)-2.*uprev(2:X-1)+uprev(1:X-2))+0.5*(dt^2).*f(2,2:X-1);

hist(1,:)=uprev;
hist(2,:)=u;

for n=3:T       
    unext(2:X-1)=-uprev(2:X-1)+2.*u(2:X-1)+(C^2).*(u(3:X)-2.*u(2:X-1)+u(1:X-2))+(dt^2).*f(n,2:X-1);
    unext(1)=0;
    unext(X)=0;
    uprev=u;
    u=unext;
    
    hist(n,:)=unext; 
end
                                                                           
fnct=hist;
end