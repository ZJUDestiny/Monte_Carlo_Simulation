function energy_res_back=energy_distributionR(r,del,lamda,energy_res_back,z,depth,hengzuobiao)
%������ӵ������������������Ϊ��ǰ���ӵ�����r dE/dS, lamda,
%���е�������������energy_res��zΪ���������꣬depthΪͳ�Ƶ������ֵ

%ͳ�Ƴ��ִ�����r��Ӧ��Ϊ�����꣬���������ֵ����
    for i=1:length(r)
        temp=floor(r(i))+1;
        if temp>=1 && z(i)<depth && temp<hengzuobiao
            energy_res_back(floor(r(i))+1)=energy_res_back(floor(r(i))+1)+del(i)*lamda(i)/pi/(floor(r(i))+1+floor(r(i)))/depth;
        end
    end
end