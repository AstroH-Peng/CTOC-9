function rv =  E2C(ele,mu)
%�������ת��Ϊ�ѿ�������
d2r = pi/180;

%��ȡ�������
a = ele(1);
e = ele(2);
i = ele(3)*d2r;
Omega = ele(4)*d2r;
omega = ele(5)*d2r;
f = ele(6)*d2r;

%��ƫ�����E��ƽ�����M

E = 2*atan2(sin(f/2)*sqrt(1-e),cos(f/2)*sqrt(1+e));
E = mod(E,2*pi);

M = E-e*sin(E);
M = mod(M,2*pi);

%��P��Q
Q = [-cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(i);
    -sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(i);
    cos(omega)*sin(i)];
P = [cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(i);
    sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(i);
    sin(omega)*sin(i)];



%�� r\v

r = a*(cos(E)-e)*P  + a*sqrt(1-e^2)*sin(E)*Q;
v = sqrt(mu*a)/norm(r)*(-sin(E)*P  +  sqrt(1-e^2)*cos(E)*Q);
rv = [r;v]';

