function rv =  E2C(ele,mu)
%轨道根数转化为笛卡尔坐标
d2r = pi/180;

%读取轨道根数
a = ele(1);
e = ele(2);
i = ele(3)*d2r;
Omega = ele(4)*d2r;
omega = ele(5)*d2r;
f = ele(6)*d2r;

%求偏近点角E和平近点角M

E = 2*atan2(sin(f/2)*sqrt(1-e),cos(f/2)*sqrt(1+e));
E = mod(E,2*pi);

M = E-e*sin(E);
M = mod(M,2*pi);

%求P和Q
Q = [-cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(i);
    -sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(i);
    cos(omega)*sin(i)];
P = [cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(i);
    sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(i);
    sin(omega)*sin(i)];



%求 r\v

r = a*(cos(E)-e)*P  + a*sqrt(1-e^2)*sin(E)*Q;
v = sqrt(mu*a)/norm(r)*(-sin(E)*P  +  sqrt(1-e^2)*cos(E)*Q);
rv = [r;v]';

