e2=1.44; % e^2=1.44 MeV fm
Zp=2; 
Zt=79; Nt=118; At=Zt+Nt;
Rt=1.2*At^(1/3); % Radius of target in fm
KEp=5.0; % KE of projectile in MeV
x0=-100;
hold on; axis equal; axis([-100 100 0 200]);
t=linspace(0,2*pi,100); [x,y]=pol2cart(t,Rt); fill(x,y,'g');
for b=0:5:200;
    k=Zp*Zt*e2./(2*KEp*b^2);
    t0=pi/2+atan(b*k);
    u0=-k/cos(t0);
    t=linspace(atan2(b,x0),2*t0-pi,100);
    u=u0*cos(t-t0)-k;
    [x,y]=pol2cart(t,1./u);
    plot(x,y)
end