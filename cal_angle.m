function [theta,fai]=cal_angle(beta)
    %ɢ���theta����λ��fai�ļ��㺯��
    R=rand();
    costheta=1-2*beta*R/(1+beta-R);
    %��ֹcosֵС��-1���������
    if costheta<-1
        costheta=-1;
    end
    theta=acos(costheta); %������
    R=rand();
    fai=2*pi*R;
end