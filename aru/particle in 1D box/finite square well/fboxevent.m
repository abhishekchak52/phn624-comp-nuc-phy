function [zdot,isterminal,dircn] = fboxevent(l,z,flag)
global e;
if nargin<3||isempty(flag)   
    c=3*10^8*10^(-15);
    h=197.327/c;
    m=938.272/c^2;
    if l<=-10 || l>=10
       V=0;
    else
       V=25;
    end
zdot=[z(2); (2*m*(e-V)/h^2)*z(1)];
else
    if flag=='events'
        zdot=z(1);
        isterminal=0;
        dircn=0;
    end
end