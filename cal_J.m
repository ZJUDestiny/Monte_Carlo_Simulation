function Jres=cal_J(Z)
%���еļ�������ܵĺ������˴���������
    Jres=zeros(size(Z,1),1);
    for i=1:size(Z,1)
        if Z(i)<=12
            Jres(i)=11.5*Z(i)*10^-3;
        else
            Jres(i)=(9.76*Z(i)+58.8/Z(i)^0.19)*10^(-3);
        end
    end
end