function [zdot,isterminal,dircn] = fboxeventws(l,z,flag)
global e;
if nargin<3||isempty(flag)   
    c=3*10^8*10^(-15);
    h=197.327/c;
    m=938.272/c^2;
    if l<=-10
        V=0;
    elseif l>=-10 && l<=0
        V=25/(1+exp((-9-l)/0.5));
    elseif l>=0 && l<=10
        V=25/(1+exp((l-9)/0.5));
    else
        V=0;
    end
zdot=[z(2); (2*m*(e-V)/h^2)*z(1)];
else
    if flag=='events'
        zdot=z(1);
        isterminal=0;
        dircn=0;
    end
end

% TO PLOT POTENTIAL
% l=linspace(-12,12,100);
% for i=1:100
%     V(i)=0;
%     if l(i)<=-10
%         V(i)=0;
%     elseif l(i)>=-10 && l(i)<=0
%         V(i)=-50/(1+exp((-5-l(i))/0.5));
%     elseif l(i)>=0 && l(i)<=10
%         V(i)=-50/(1+exp((l(i)-5)/0.5));
%     else
%         V(i)=0;
%     end
% end
% plot(l,V)

% l1=linspace(-10,0,100)
% r=length(l1);
% for i=1:r
%  V1(i)=-50/(1+exp((-5-l1(i))/0.5));
% end
% plot(l1,V1,'k');
% hold on
% l2=linspace(0,10,100)
% s=length(l2);
% for i=1:s
%  V2(i)=-50/(1+exp((l2(i)-5)/0.5));
% end
% plot(l2,V2,'k');