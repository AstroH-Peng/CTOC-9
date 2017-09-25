function [out]=modelJ2(t,x)
%J2摄动轨道动力学模型
global mu;
global J2;
global Re;


out=zeros(6,1);                       

r=sqrt(x(1)^2+x(2)^2+x(3)^2);                               %探测器距地球距离

out(1)=x(4);
out(2)=x(5);
out(3)=x(6);
out(4)=-mu*x(1)/r^3*(1+3/2*J2*(Re/r)^2*(1-5*x(3)^2/r^2));
out(5)=-mu*x(2)/r^3*(1+3/2*J2*(Re/r)^2*(1-5*x(3)^2/r^2));
out(6)=-mu*x(3)/r^3*(1+3/2*J2*(Re/r)^2*(3-5*x(3)^2/r^2));

end