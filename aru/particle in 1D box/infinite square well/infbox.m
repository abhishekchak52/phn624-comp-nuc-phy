hold on
global e;
lspan=[-10 10];
e=0; einc=-1.0;
z0=[0; 0.01];
options=odeset('events','on');
[l,z,le,ze,ie]=ode45('infbox_ode',lspan,z0,options);
ie1=length(ie);
E1=e;
n=0;
Ef=[];
while n < 5
    e=e+einc;
    [l,z,le,ze,ie]=ode45('infbox_ode',lspan,z0,options);
    ie2=length(ie);
    E2=e;    
    if ie2==ie1+1
        for s=1:400
            e=(E1+E2)/2;
            [l,z,le,ze,ie]=ode45('infbox_ode',lspan,z0,options);
            ie3=length(ie);
            E=e;
            if ie3==ie2
                E2=E;
            else
                E1=E;
            end
            if(abs(E1-E2)<1e-5)
                break;
            end
        end
        Ef=[Ef E];
        plot(l,z(:,1));
        n=n+1;
    end
    E1=E2;ie1=ie2;
end
legend(cellstr(num2str(Ef','E = %10.4f'))');
xlabel('x');ylabel('\psi (x)');title('Infinite Well');
