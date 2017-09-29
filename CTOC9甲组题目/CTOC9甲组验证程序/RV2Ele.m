function Ele = RV2Ele(RV,mu)
%==========================================================================
%   ���ֳ�   lilincheng2012@163.com      2017-09-10
%==========================================================================
% ���ܣ������ǵ�λ���ٶ�ת��Ϊ�������,��������������
%
% ���룺
%      mu             ������������ km2/s2
%      rv             ���ǵ�λ���ٶȣ���λ��km, km/s
% �����
%      Ele            ���ǵĹ��������[a(km) e i(rad) RAAN(rad) omg(rad) f(rad)]
%
% ע�⣺               λ���ٶȺ͹������������������ʽ����
%
% �ο����ף�           Howard D. Curtis, Orbit Mechanics[M]
%==========================================================================

R  = RV(:,1:3);
V  = RV(:,4:6);
No = size(RV,1);

r = sqrt(RV(:,1).^2 + RV(:,2).^2 + RV(:,3).^2);
v = sqrt(RV(:,4).^2 + RV(:,5).^2 + RV(:,6).^2);
vr = dot(R',V') ./ r';
H  = cross(R,V);
h = sqrt(H(:,1).^2 + H(:,2).^2 + H(:,3).^2);
N = cross(repmat([0 0 1],No,1),H);
n = sqrt(N(:,1).^2 + N(:,2).^2 + N(:,3).^2);

% ��������incl
incl = acos(H(:,3) ./ h);

% ���������ྭRAAN
for ii =1 : size(RV,1)
    if n(ii) ~= 0
        RAAN(ii,1) = acos(N(ii,1) / n(ii));
        if N(ii,2) < 0
            RAAN(ii,1) = 2 * pi - RAAN(ii,1);
        end
    else
        RAAN(ii,1) = 0;
    end
end

% ���ƫ����e
A = repmat((v.^2  - mu ./ r),1,3);
B = repmat((r .* vr'),1,3);
E = 1 / mu * (A .* R - B .* V);
for ii =1 : size(RV,1)
    e(ii,1) = norm(E(ii,:));
end

% �����ص����omg
for ii =1 : size(RV,1)
    if n(ii) ~= 0
        if e(ii,1) > eps
            omg(ii,1) = acos( dot( N(ii,:),E(ii,:) ) / n(ii) / e(ii,1) );
            omg(ii,1) = real(omg(ii,1));
            if E(ii,3) < 0
                omg(ii,1) = 2 * pi - omg(ii,1);
            end
        else
            omg(ii,1) = 0;
        end
    else
        omg(ii,1) = 0;
    end
    
% ���������f
    if e(ii,1) > eps
        f(ii,1) = acos(dot(E(ii,:),R(ii,:)) / e(ii,1) / r(ii,1));
        f(ii,1) = real(f(ii,1));
        if vr(ii) < 0
            f(ii,1) = 2 * pi - f(ii,1);
        end
    else
        cp = cross(N(ii,:),R(ii,:));
        if cp(3) >= 0
            f(ii,1) = acos( dot( N(ii,:),R(ii,:) ) / n(ii) / r(ii,1) );
        else
            f(ii,1) = 2 * pi - acos( dot( N(ii,:),R(ii,:) ) / n(ii) / r(ii,1) );
        end
    end
end

% ���볤��a
a   = h.^2 ./ mu ./ (1 - e.^2);

% ��������ྭ�����ص���Ǻ�������ͶӰ��[0,2*pi]֮��
RAAN = mod( RAAN , 2 * pi );
omg  = mod( omg  , 2 * pi );
f    = mod( f    , 2 * pi );

% ����������
Ele = [a e incl RAAN omg f];



