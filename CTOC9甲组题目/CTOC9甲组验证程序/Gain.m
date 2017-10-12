function GainComp = Gain(TaskNo, TN, RVN1, RVN2, RVN3, lamda, El, Az) 
%==========================================================================
%   ����   zenghaohz@163.com      2017-09-10
%==========================================================================
% ���ܣ���ǰ�������ж�GEO���ǵļ������
% ���룺 
%      TaskNo         ������
%      TN             ����ʱ�����У�Nά����������λ��s
%      RVN1           �����ڶ�Ӧ��TN�ļ������1��״̬���У�N*6ά���󣬵�λ��km��km/s
%      RVN2           �����ڶ�Ӧ��TN�ļ������2��״̬���У�N*6ά���󣬵�λ��km��km/s
%      RVN3           �����ڶ�Ӧ��TN�ļ������3��״̬���У�N*6ά���󣬵�λ��km��km/s
%      lamda          GEO���Ƕ��㾭�ȣ�deg
%      El             ����ָ�����ǣ�deg
%      Az             ����ָ��λ�ǣ�deg
% �����
%      GainComp       ��ǰ�������ж�GEO���ǵļ�����棬��λ��km^2
%==========================================================================

RN = [RVN1(:,1:3)'; RVN2(:,1:3)'; RVN3(:,1:3)'];
GainN = zeros(length(TN),1);

% ����
rGEO = 42164.2; % GEO����뾶����λ��km
Thetag0 = 100.1218 * pi / 180; % ��ʼʱ��(MJD = 58849)�ĸ������κ���ʱ�ǣ���λ��rad
Dthetag = 360 / 0.99727 / 86400 * pi / 180; % �������κ���ʱ�ǵı仯�ʣ���λ��rad/sec

% ��λ���㣺deg -> rad
lamda = lamda * pi / 180;
El = El * pi / 180;
Az = Az * pi / 180;

for jj = 1:length(TN)
    
    % GEO������J2000����ϵ�µ�λ��ʸ��
    Thetag = Thetag0 + Dthetag * TN(jj); %rad
    Mie = [ cos(-Thetag)  sin(-Thetag) 0;
           -sin(-Thetag)  cos(-Thetag) 0;
            0            0           1];
    rGEO_ECF = [rGEO * cos(lamda); rGEO * sin(lamda); 0];
    rGEO_J2000 = Mie * rGEO_ECF;

    % GEO���ǲ���ʸ����J2000����ϵ�µı�ʾ
    g = [cos(El) * cos(Az); cos(El) * sin(Az); sin(El)]; % % ����ʸ�����������ϵ��
    Mx90 = [1,  0,               0;
            0,  cos(90*pi/180),  sin(90*pi/180);
            0, -sin(90*pi/180),  cos(90*pi/180)];
    Mz90 = [ cos(-90*pi/180), sin(-90*pi/180), 0;
            -sin(-90*pi/180), cos(-90*pi/180), 0;
             0,               0,               1];
    MzTL = [ cos(- Thetag - lamda), sin(- Thetag - lamda), 0;
            -sin(- Thetag - lamda), cos(- Thetag - lamda), 0;
             0,                     0,                     1];
    M     = MzTL * Mz90 * Mx90;
    Beam = M * g; % ����ʸ����J2000����ϵ��
    
    for kk = 1 : 3
        
        % ������������GEO��λ��ʸ��
        Ri = RN(1+(kk-1)*3 : 3+(kk-1)*3,jj);
        Rgi = Ri - rGEO_J2000;

        % ����������Ƿ��ڲ�����
        Angle1 = acos( dot(Rgi, Beam) / (norm(Rgi) * norm(Beam)));
%         if Angle1 * 180 / pi > 0.5 
%             fprintf('�� %d ���������� �� %d �ż�����ǲ��ڲ���Բ׶��\n',TaskNo, kk);
%         end
        
        % �����������GEO�����Ƿ��ڵ���ͬ�� 
        Angle2 = acos( dot(Ri, rGEO_J2000) / (norm(Ri) * norm(rGEO_J2000)));
%         if Angle2 * 180 / pi  > 90
%             fprintf('�� %d ���������� �� %d �ż��������GEO���ǲ���ͬ��\n',TaskNo, kk);
%         end        
    end
    
    % ����������
    R12 = RN(1:3,jj) - RN(4:6,jj);
    R23 = RN(4:6,jj) - RN(7:9,jj);
    R13 = RN(1:3,jj) - RN(7:9,jj);
    L = (norm(R12) + norm(R13) + norm(R23))/2;
    S = sqrt( L * (L - norm(R12)) * ( L- norm(R13)) * (L - norm(R23)) );
    temp = cross(R12, R13);
    gk = abs(dot(Beam, temp/norm(temp))); % ������������ƽ��нǵ�����
    GainN(jj) = S * gk;     
end

% ��ǰ������ļ��������
GainComp = sum(GainN);

end

