function energy_res=energy_distribution(z,del,lamda,zmax,energy_res)
%������ӵ������������������Ϊ��ǰ���ӵ�����z dE/dS, lamda, ���������zmax��
%���е�������������energy_res��zΪ���������꣬depthΪͳ�Ƶ������ֵ

%ͳ�Ƴ��ִ�����z��Ӧ��Ϊ�����꣬���������ֵ����
    for i=1:length(z)
        temp=floor(z(i))+1;
        if temp<=zmax && temp>1
            energy_res(floor(z(i))+1)=energy_res(floor(z(i))+1)+del(i)*lamda(i);
        end
    end
end