clear
close all
tic

%����һ��cell�������洢���еĲ���
layer={};
%������һ��cell�������洢��ά���������Լ�������Ϣ
electron={};
hengzuobiao=150;
%�û�ѡ�����뷽ʽ
type=input('��ѡ��������루0�������ļ����루1��: ');
if type==0  
    %��������
    %��Ҫ����ĵ�����Ŀ
    count_ele=input('���������ĵ�����Ŀ��');
    %��������ȷ��������ֲ�ͼ��ķ�Χ��z�����ֵ��
    zmax=input('�������������ȷ��������ֲ�ͼ��ķ�Χ��z�����ֵ����');
    E=input('���������������ֵ');
    n=input('���������: ');
    for i=1:n
        %Ci
        C_temp=input(['�������',num2str(i),'����ϵ�Ci: ']);   
        layer{i,1}=C_temp;
        %Zi
        Z_temp=input(['�������',num2str(i),'����ϵ�Zi: ']);   
        layer{i,2}=Z_temp;
        %Ai
        A_temp=input(['�������',num2str(i),'����ϵ�Ai: ']);
        layer{i,3}=A_temp;
        %rho
        rho_temp=input(['�������',num2str(i),'����ϵ��ܶ�: ']);
        layer{i,4}=rho_temp;
        %depth
        zmax_temp=input(['�������',num2str(i),'����ϵĺ��: ']);
        layer{i,5}=zmax_temp;
        %name
        name=input(['�������',num2str(i),'����ϵ�����: '],'s');
        layer{i,6}=name;
    end
else
    %�ļ�����
    disp('�ļ��ڵ�һ��Ϊ���������Ŀ���ڶ���Ϊ��ͼʱZ����ȣ�������Ϊ��������E��֮��ÿһ��Ϊ������˳��ΪCi Zi Ai rho depth');
    %���ַ�����ʽ��ȡ����������ļ���
    dataname=input('�������ļ���: ','s'); 
    file=fopen(dataname,'r');
    %��Ҫ����ĵ�����Ŀ
    count_ele=str2double(fgetl(file)); 
    %��������ȷ��������ֲ�ͼ��ķ�Χ��z�����ֵ��
    zmax=str2double(fgetl(file)); 
    %��õ�������E
    E=str2double(fgetl(file)); 
    %��ʼ��n
    n=0; 
    %���Ƿ��ļ���βΪ�ж��Ƿ��ȡ����
    while ~feof(file)
        n=n+1;
        line_now=split(fgetl(file));    %��ÿһ��������Ƭ
        layer{n,1}=eval(line_now{1});   %Ci
        layer{n,2}=eval(line_now{2});   %Zi
        layer{n,3}=eval(line_now{3});   %Ai
        layer{n,4}=eval(line_now{4});   %rho
        layer{n,5}=eval(line_now{5});   %depth       
        layer{n,6}=line_now{6};   %name
    end
    fclose(file);
end
%�ļ�������"depth_E_count_name1_name2"
filepath=[num2str(layer{1,5}),'_',num2str(E),'_',num2str(count_ele),'_',layer{1,6},'_',layer{2,6}];

%�ж��ļ����Ƿ���ڣ�������������½����ļ���
if ~exist(['.\','result\',filepath],'dir')
    mkdir(['.\','result\',filepath]);
else
    str_temp=input('Ĭ���ļ����Ѵ��ڣ��������ļ�������ĩβ���ֱ�ʶ�����������룬�򸲸�ԭ�ļ���');
    filepath=[filepath,num2str(str_temp)];
    mkdir(['.\','result\',filepath]);
end

%��ÿһ��ĺ�ȴ���һ����Ϊdepth��������
depth=zeros(1,n);
for n_layer=1:n
    depth(n_layer)=layer{n_layer,5};
end

%��ʼ�������������󡢱�ɢ�����������������͸�����������������
energy_res=zeros(1,floor(zmax)+1);
energy_res_back=zeros(1,hengzuobiao);
energy_res_forward=zeros(1,hengzuobiao);

%����ˮƽ�ߵĺ����귶Χ
min_y=0;max_y=0;
min_x=0;max_x=0;
minx=[];maxx=[];miny=[];maxy=[];minz=[];maxz=[];
%��ɢ����ӵ���Ŀ
count_back=0;

%��ÿһ�����ӽ���ѭ������
for jjj=1:count_ele
    %ÿ��100���������һ�Σ����û�Ԥ��ʱ��
    if mod(jjj,100)==0
        disp([num2str(jjj/100),' / ',num2str(count_ele/100)]);
    end
    
    %����code_1���㲢���Ƶ��ӵ���ϸ�켣��Ϣ
    [energy_res_back,energy_res_forward,energy_res,miny_now,maxy_now,electron_temp,coor_temp,minx_now,maxx_now,count_back]=code_1(energy_res,zmax,E,layer,depth,energy_res_back,energy_res_forward,count_back,hengzuobiao);
    %������ֵ
    if miny_now<min_y
        min_y=miny_now;
    end
    if maxy_now>max_y
        max_y=maxy_now;
    end
    
    if minx_now<min_x
        min_x=minx_now;
    end
    if maxx_now>max_x
        max_x=maxx_now;
    end
    
    %�������ֲ���άͼ���������ֵ���Լ�ÿ�����ӵ���������
    minx=[minx,coor_temp(1)];maxx=[maxx,coor_temp(2)];
    miny=[miny,coor_temp(3)];maxy=[maxy,coor_temp(4)];
    minz=[minz,coor_temp(5)];maxz=[maxz,coor_temp(6)];
    electron{jjj}=electron_temp;
end
%��άͼ�����ֵ
Min_x=min(minx);Max_x=max(maxx);
Min_y=min(miny);Max_y=max(maxy);
Min_z=min(minz);Max_z=max(maxz);

%��ֵ�и����Լ�С������Ϊ������࣬���ҳ�ʼ����ά����
deltax=ceil(Max_x-Min_x+1);
deltay=ceil(Max_y-Min_y+1);
deltaz=ceil(Max_z-Min_z+1);
energy_distribution3=zeros(deltay,deltaz,1);

%��ÿһ�����ӽ��в���
for k=1:count_ele
    %��ȡ�õ�������  
    electron_temp=electron{k};
    for kk=1:length(electron_temp{1})
        %����ά�����Ӧ��λ�õ�������õ��������������
        temp=energy_distribution3(ceil(electron_temp{2}(kk)+1-Min_y),ceil(electron_temp{3}(kk)+1-Min_z),1);
        temp=temp+electron_temp{4};
        energy_distribution3(ceil(electron_temp{2}(kk)+1-Min_y),ceil(electron_temp{3}(kk)+1-Min_z),1)=temp;   
    end

end
%��Ϊ��ֵ�������ͼ��
energy_distribution3=-energy_distribution3;

%����ֵ������ֵ256���Ի��ƻҶ�ͼ
energy_distribution3=256/max(max(energy_distribution3))*energy_distribution3;

%��ÿ��С��Ԫ�������ۼӣ�delta_y��delta_zΪС��Ԫ�ı߳�
delta_y=40;
delta_z=40;
%С��Ԫ�����Ŀ
count_y=floor(deltay/delta_y); 
count_z=floor(deltaz/delta_z);
%��ÿһ��С��Ԫ����в��������ڲ�Ԫ�صĺ���Ϊ�ڲ�����Ԫ���µ�ֵ
for kx=1:count_y
    for kz=1:count_z
        energy_distribution3((kx-1)*delta_y+1:(kx)*delta_y,(kz-1)*delta_z+1:kz*delta_z)=ones(delta_y,delta_z).*sum(sum(energy_distribution3((kx-1)*delta_y+1:(kx)*delta_y,(kz-1)*delta_z+1:kz*delta_z))); 
    end
end

%��ȥ��Ե�����������Ԫ��ֵ
energy_distribution3=energy_distribution3(1:delta_y*count_y,1:delta_z*count_z);
%�����ֵ��Ϊ256���Ի��ƻҶ�ͼ
energy_distribution3=256/max(max(energy_distribution3))*energy_distribution3;

%��ת��ʹ�ñ���Ϊ��ɫ
energy_distribution3=-1.*energy_distribution3+256;
%energy_distribution3_back=-1.*energy_distribution3_back+256;
%����洢��ʱ�򣬺������귴�ˣ��˴�����ת��
energy_distribution3=energy_distribution3';

%��ͼ���������������Ķ�ά�Ҷ�ͼ
figure(2);imagesc(energy_distribution3);colormap(gray(256));title('energy distribution');

%ָ�����»��ƣ��ʽ�zȡ��
energy_res=-energy_res;
energy_res_back=-energy_res_back;
energy_res_forward=-energy_res_forward;

%�ڵ��ӹ켣��ͼ���ϻ��Ʋ��ϵķֽ���
for k=1:n-1
    if k==1
    %����ǵ�һ���ֽ��ߣ��ڴ˴��ֽ���
        figure(5);line([min_x-200,max_x+200],[-(layer{k,5}-1),-(layer{k,5}+1)]);
        figure(1);[X,Y,z] = meshgrid(min_x-200:30:max_x+200,min_y-200:30:max_y+200,-(layer{k,5}));
    else
        figure(5);line([min_x-200,max_x+200],[-(layer{k,5}+layer{k-1,5}-1),-(layer{k,5}+layer{k-1,5}+1)]);
        figure(1);[X,Y,z] = meshgrid(min_x-200:30:max_x+200,min_y-200:30:max_y+200,-(layer{k,5}+layer{k-1,5}));
    end
end
mesh(X,Y,z)

%������Z���Լ�R����Ƶ�ͼ��ĺ�������
energy_his=zeros(1,length(energy_res));
energy_his_back=zeros(1,length(energy_res_back));
energy_his_forward=zeros(1,length(energy_res_forward));

%ָ������ȷ������������ĺ�������Ϊ50nm����R�������������ĺ�������Ϊ20nm
delta_his=1;
delta_his_R=1;
%��ÿһ��delta_his�����������ۼӣ���Ϊ�ⲿ�ֵ�ֵ
for t=1:floor(zmax/delta_his)
    energy_his((t-1)*delta_his+1)=sum(energy_res((t-1)*delta_his+1:t*delta_his))/(count_ele*E);
end
for t=1:hengzuobiao/delta_his_R
    energy_his_back((t-1)*delta_his_R+1)=sum(energy_res_back((t-1)*delta_his_R+1:t*delta_his_R));
    energy_his_forward((t-1)*delta_his_R+1)=sum(energy_res_forward((t-1)*delta_his_R+1:t*delta_his_R));
end

for i=1:length(energy_his)
    figure(6);
    tempzong=linspace(0,energy_his(i),2);
    plot(i.*ones(length(tempzong)),tempzong,'color','b');hold on;
end
%������ȷ������������ͼ
figure(3);plot([1:zmax+1],energy_his);title('energy distribution histogram');axis([-50,zmax+100,0,1.2*max(energy_his)]);xlabel('z/nm');
%����ɢ����ӡ�͸�������R�������������ͼ
figure(4);subplot(131);plot([1:hengzuobiao],energy_his_back);title('back');axis([-10,hengzuobiao,0,1.2*max(energy_his_back)]);xlabel('R/nm');
figure(4);subplot(132);plot([1:hengzuobiao],energy_his_forward);title('forward');axis([-10,hengzuobiao,0,1.2*max(energy_his_forward)]);xlabel('R/nm');
figure(4);subplot(133);plot([1:hengzuobiao],energy_his_back+energy_his_forward);title('all');axis([-10,hengzuobiao,0,1.2*max(energy_his_back+energy_his_forward)]);xlabel('R/nm');

%�½�һ����¼��ɢ��ϵ�����ı��ļ�������ɢ��ϵ����¼������
backfile=fopen(['.\','result\',filepath,'\backscattered.txt'],'w');
fwrite(backfile,['Backscattering cofficient: ',num2str(count_back/count_ele)]);
fclose(backfile);
disp(['Backscattering cofficient: ',num2str(count_back/count_ele)]);

%������ͼ�񶼱��浽����������ļ���֮�У���ʽΪpng��ʽ��fig��ʽ����
saveas(1,['.\','result\',filepath,'\trajectory3D.png']);
saveas(1,['.\','result\',filepath,'\trajectory3D.fig']);
saveas(2,['.\','result\',filepath,'\energy2D.png']);
saveas(2,['.\','result\',filepath,'\energy2D.fig']);
saveas(3,['.\','result\',filepath,'\energy1D.png']);
saveas(3,['.\','result\',filepath,'\energy1D.fig']);
saveas(5,['.\','result\',filepath,'\trajectory2D.png']);
saveas(5,['.\','result\',filepath,'\trajectory2D.fig']);
saveas(4,['.\','result\',filepath,'\energy_radius.png']);
saveas(4,['.\','result\',filepath,'\energy_radius.fig']);
saveas(6,['.\','result\',filepath,'\energy_z.png']);
saveas(6,['.\','result\',filepath,'\energy_z.fig']);
toc