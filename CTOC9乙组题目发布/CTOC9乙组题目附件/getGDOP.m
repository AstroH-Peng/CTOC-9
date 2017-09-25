% =========================================================================
% ����GDOP
%
% ����:
% satMtx    [n x 3]    �ɼ�������ecfϵ�µ�λ������
% obsEcf    [1 x 3]    ������ecfϵ�µ�λ������
%
% ���:
% gdop      [1 x 1]     ���ξ�������
%
% ���ԣ�
% satMtx = [15524.471175 , -16649.826222 , 13512.2723887;
%     -2304.058534 , -23287.906465 , 11917.038105;
%     16680.243357 , -3069.625561 , 20378.551047;
%     -14799.931395 , -21425.35824 , 6069.947224];
% obsEcf = [-730.000 , -5440.000 , 3230.000];
% gdop_(satMtx , obsEcf)
%
% �ο����ף�
% https://en.wikipedia.org/wiki/Dilution_of_precision_(navigation)
%
% Copyright �ų�(Chen Zhang)
% �廪��ѧ�����(Tsinghua Space Center)
% chenzhang86@tsinghua.edu.cn
% 2017/09/01
% -------------------------------------------------------------------------
function gdop = getGDOP(satMtx , obsEcf)
seeNum = size(satMtx , 1); % �ɼ�������
rTemp = satMtx - repmat(obsEcf , seeNum , 1); % ����վָ�����ǵ�����
rNorm = rTemp ./ repmat(sqrt(rTemp(: , 1).^2 + rTemp(: , 2).^2 + rTemp(: , 3).^2) , 1 , 3); % ����վָ�����ǵĵ�λ����
H = [rNorm , ones(seeNum , 1)]; % ����H����
gdop = sqrt( trace(pinv(H' * H)) ); % ����gdop
end
