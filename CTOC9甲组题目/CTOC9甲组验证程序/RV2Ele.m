function Ele = RV2Ele(RV,mu)
%==========================================================================
%   李林澄   lilincheng2012@163.com      2017-09-10
%==========================================================================
% 功能：将卫星的位置速度转化为轨道根数,允许多组输入输出
%
% 输入：
%      mu             地球引力常数 km2/s2
%      rv             卫星的位置速度，单位：km, km/s
% 输出：
%      Ele            卫星的轨道根数，[a(km) e i(rad) RAAN(rad) omg(rad) f(rad)]
%
% 注意：               位置速度和轨道根数均以行向量形式给出
%
% 参考文献：           Howard D. Curtis, Orbit Mechanics[M]
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

% 求解轨道倾角incl
incl = acos(H(:,3) ./ h);

% 求解升交点赤经RAAN
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

% 求解偏心率e
A = repmat((v.^2  - mu ./ r),1,3);
B = repmat((r .* vr'),1,3);
E = 1 / mu * (A .* R - B .* V);
for ii =1 : size(RV,1)
    e(ii,1) = norm(E(ii,:));
end

% 求解近地点幅角omg
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
    
% 求解真近点角f
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

% 求解半长轴a
a   = h.^2 ./ mu ./ (1 - e.^2);

% 将升交点赤经，近地点幅角和真近点角投影在[0,2*pi]之间
RAAN = mod( RAAN , 2 * pi );
omg  = mod( omg  , 2 * pi );
f    = mod( f    , 2 * pi );

% 输出轨道根数
Ele = [a e incl RAAN omg f];



