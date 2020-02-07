set parametric
set urange [0:2*pi]
set vrange [-pi/2:pi/2]
r=2
fx(v,u)=r*cos(v)*cos(u)
fy(v,u)=r*cos(v)*sin(u)
fz(v)=r*sin(v)
set view equal xyz
set isosamples 100,100; set samples 100,100
splot fx(v,u),fy(v,u),fz(v) w pm3d
pause -1
