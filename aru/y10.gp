set parametric
set urange [0:2*pi]
set vrange [-pi/2:pi/2]
fx(v,u)=cos(v)*cos(v)*cos(u)
fy(v,u)=cos(v)*cos(v)*sin(u)
fz(v)=cos(v)*sin(v)
set view equal xyz
set isosamples 300,300; set samples 300,300
splot fx(v,u),fy(v,u),fz(v) w pm3d
pause -1
