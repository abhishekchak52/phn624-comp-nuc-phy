hold on
global e;
lspan=[-12 12];
Ei=0:25;
Ei=25-Ei;
r=length(Ei)-1;
z0=[0; 0.01];
for i=1:r
    e=Ei(i);
    options=odeset('events','on');
    [l,z,le,ze,ie]=ode45('fboxevent',lspan,z0,options);
    ie1=length(ie);
    E1=e;
    j=i+1;
    e=Ei(j);
    options=odeset('events','on');
    [l,z,le,ze,ie]=ode45('fboxevent',lspan,z0,options);
    ie2=length(ie);
    E2=e;    
    if ie2==ie1+1
        for s=1:50
            e=(E1+E2)/2;
            options=odeset('events','on');
            [l,z,le,ze,ie]=ode45('fboxevent',lspan,z0,options);
            ie3=length(ie);
            E=e;
            if ie3==ie2
                E2=E;
            else
                E1=E;
            end
        end
        Ef=25-E
        x=z(:,1);
        plot(l,x);
    else
    end
end

q=linspace(-12,12,1000);
plot(q,0,'k');
w=linspace(-0.2,0.32,1000);
plot(-10,w,'k',10,w,'k');
axis([-12 12 -0.2 0.32]);