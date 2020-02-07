global e2 mp Rt qQ;
e2=1.44;  % e^2=1.44 MeV fm
mn=939.5654133; % Mass of neutron in MeV/c^2
mz=938.2720813; % Mass of proton in MeV/c^2
Zp=2; Np=2; Ap=Zp+Np; mp=Zp*mz+Np*mn;
Zt=79; Nt=118; At=Zt+Nt; mt=Zt*mz+Nt*mn;
qQ=Zp*Zt;
Rt=10.2*At^(1/3); % Radius of target in fm
KEp=5.0; % KE of projectile in MeV
vp=sqrt(2*KEp/mp);
hold on; axis equal; axis([-100 200 -200 200]);
t=linspace(0,2*pi,100); [x,y]=pol2cart(t,Rt); fill(x,y,'g');
for b=-200:5:200;
    tspan=linspace(0,10000,100);  % Time in units of fm/c
    [t,x]=ode45(@rscatode,tspan,[-100; vp; b; 0]);
    plot(x(:,1),x(:,3))
end