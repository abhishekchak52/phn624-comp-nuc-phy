function [dz] = rscatode(t,z)
global e2 mp Rt qQ;
r=sqrt(z(1).^2+z(3).^2);
r=max(r,Rt);
dz=[z(2); qQ*e2.*z(1)/(mp*r.^3); z(4); qQ*e2.*z(3)/(mp*r.^3)];
