function [energy_res_back,energy_res_forward,energy_res,miny,maxy,electron_temp,coor_temp,minx,maxx,count_back]=code_1(energy_res,zmax,E,layer,depth_all,energy_res_back,energy_res_forward,count_back,hengzuobiao)
%������������е������ֲ�����energy_res��z�����ֵzmax������������ʼ����E���������в�����layer�����в�����depth_all
%��ɢ����ӵ�������������energy_res_back��͸����ӵ�������������energy_res_forward����ɢ����ӵ���Ŀcount_back

%������������ϴ˵��ӵ������ֲ�����energy_res_back, energy_res_forward, energy_res��
%�˵��ӵ��������Сֵminx miny�����ֵmaxx maxy����ǰ���ӵ���������electron_temp����������electron_temp��
%��ɢ�������Ŀcount_back

    %�����Ƿ�ɢ��
    is_back=0;
    
    %��һ�ε�depth�������ۼ�depth�������Ƿ�Խ���߽�
    depth_last=0;
    %��ʼ�������ֲ����������
    %ÿһ�����ӵ�x,y,z����
    x_x=[];
    y_x=[];
    z_x=[];
    %ÿһ�����Ӷ�Ӧ��z��dE/dS
    del_energy=[];
    %ÿһ�����Ӷ�Ӧ��z��lamda
    lamda_energy=[];
    
    %��ʼ����������λ�õľ���
    finalx=[];
    finaly=[];
    finalz=[];

    %����٤���޳���
    NA=6.02*10^23; 

    %�����ʼ����
    %�������뾶1nm�����ѵ�λ���������
    Re=1; 
    Re=Re*10^-7; 
    %�����㵽һ�������
    w1=rand();w2=rand(); 
    x0=Re*sqrt(-log(w1))*cos(2*pi*w2);
    y0=Re*sqrt(-log(w1))*sin(2*pi*w2);
    z0=0;
    
    %����ǰ����洢����������ľ���֮�У�����10^7��Ϊ���׵�λ
    x_now=x0;y_now=y0;z_now=z0;
    finalx=[finalx,x_now*10^7];finaly=[finaly,y_now*10^7];finalz=[finalz,z_now*10^7];

    %��ȡ��һ�ֲ��ϵĲ���ֵ
    count=1;
    Ci=layer{count,1};
    Zi=layer{count,2};
    Ai=layer{count,3};
    rho=layer{count,4};
    depth=depth_last+layer{1,5};
    
    %depth_add����õ�ǰ����ȵľ���ֵ
    depth_add=zeros(1,length(depth_all)+1);
    for n_depth=1:length(depth_add)-1
        depth_add(n_depth+1)=sum(depth_all(1:n_depth));
    end
    
    %�ж��Ƿ���ײ����������ΪE>0.2keV
    while E>0.2
        %�ж��Ƿ��������һ�ֲ���֮��,�õ��������Ƿ񳬹����е�����depth֮�ͣ������ǲ��õ�ǰ���ϵ�depth
        for layer_count=1:length(depth_add)-1
            %�жϵ��ӽ�������һ��Ĳ���֮�У��������µ�һ����ϲ���
            if z_now*10^7>=depth_add(layer_count) && z_now*10^7<depth_add(layer_count+1)
                count=layer_count;
                Ci=layer{count,1};
                Zi=layer{count,2};
                Ai=layer{count,3};
                rho=layer{count,4};
                depth_last=depth;
                depth=depth_last+layer{count,5};
                break;
            end
        end   
        
        %���β�����E-keV
        beta_i=5.43*10^-3.*Zi.^(2/3)./E; 
        %�õ����ε�Rutherford����ɢ����棬��λcm2
        sigma_i=5.21*10^-21.*Zi.*(Zi+1)./E^2.*((E+511)./(E+1022))^2*4*pi./beta_i./(beta_i+1); 
        %��ԭ��ɢ�����
        ni_sigmai=Ci.*rho.*NA./Ai.*sigma_i;  
        %ƽ�����ɳ�
        lamda_ba=1/sum(ni_sigmai);
        %��ÿ��ԭ����ײ�ĸ���
        P=ni_sigmai./sum(ni_sigmai);  
        %�ж��ڵ�ǰ�����к�˭��ײ
        %�����ж����ĸ�ԭ����ײ�������R0
        R0=rand();
        %��ȡ��ǰ����ԭ��������Ŀ
        n_atom=length(Ci);
        %���С�ڵ�һ�ֵ���ײ���ʣ�i=1
        if 0<=R0 && R0<=P(1)
            i=1;
        else
        %��������жϣ����������е�ԭ�Ӹ��ʺͣ�С�ڼ�����һ��ԭ�ӵĸ��ʺͣ�����һ��ԭ����ײ��i=nn+1
            for nn=1:n_atom-1
                if R0>=sum(P(1:nn)) && R0<sum(P(1:nn+1))
                    i=nn+1;
                end
            end
        end
        %��ȡ��ײ��ԭ�ӵ�����ϵ�������Ҽ������ɳ�
        beta=beta_i(i);     
        lamda=-lamda_ba*log(rand());
        %ɢ���theta����λ��fai
        [theta,fai]=cal_angle(beta);
        
        %�����������
        %���ú�������Ji
        Ji=cal_J(Zi);
        
        %���㻯�����ƽ��ԭ��������ƽ��ԭ������ƽ��������
        %���������Ч��Z_ba A_ba J_ba�ĵ���
        Z_ba=sum(Ci.*Zi./Ai)/sum(Ci./Ai);
        A_ba=1/sum(Ci./Ai);
        J_ba=exp(sum(Ci.*Zi./Ai.*log(Ji))/sum(Ci.*Zi./Ai));
      
        %����dE/dS���Լ��������������а���������ϵ��
        del=-7.85*10^4*rho*Z_ba/A_ba/E*log(1.166*E/J_ba+0.0063*Z_ba+0.7334);
        
        %��������ת��
        %�ж����ת������Q�������ڣ���ô��ʼ��Ϊ��λ����
        if ~exist('Q')
            Q=eye(3);
        else
        %������ڣ���ô���ú�������Qn+1,n�����봫���Qn,0������˷������Qn+1,0
            Q=cal_Q(Q,theta,fai);
        end
        
        %��x,y,z���ꡢdE/dS�Լ�lamda�����������Լ��������ֲ�
        x_x=[x_x,x_now*10^7];
        y_x=[y_x,y_now*10^7];
        z_x=[z_x,z_now*10^7];
        
        %�����ϴεĵ�������
        x_last=x_now;y_last=y_now;z_last=z_now;
        %������ײ��ĵ���λ��
        [x_now,y_now,z_now]=cal_newcoor(Q,lamda,theta,fai,x_now,y_now,z_now); 
        
        
        %�����жϲ���������������Ҫ ����ȫ�²���
        for layer_count_2=1:length(depth_add)-1
            if z_now*10^7>=depth_add(layer_count_2) && z_now*10^7<depth_add(layer_count_2+1)
                %�жϲ���
                count_2=layer_count_2;
                if count_2~=count
                    %���ֿ������
                    %�µ�һ�������
                    Ci_2=layer{count_2,1};
                    Zi_2=layer{count_2,2};
                    Ai_2=layer{count_2,3};
                    rho_2=layer{count_2,4};
                    %���β�����E-keV
                    beta_i_2=5.43*10^-3.*Zi_2.^(2/3)./E; 
                    %�õ����ε�Rutherford����ɢ����棬��λcm2
                    sigma_i_2=5.21*10^-21.*Zi_2.*(Zi_2+1)./E^2.*((E+511)./(E+1022))^2*4*pi./beta_i_2./(beta_i_2+1); 
                    %��ԭ��ɢ�����
                    ni_sigmai_2=Ci_2.*rho_2.*NA./Ai_2.*sigma_i_2;  
                    %ƽ�����ɳ�
                    lamda_ba_2=1/sum(ni_sigmai_2);
                    %ɢ�䲽��
                    
                    Ji_2=cal_J(Zi_2);
                    %���㻯�����ƽ��ԭ��������ƽ��ԭ������ƽ��������
                    %���������Ч��Z_ba A_ba J_ba�ĵ���
                    Z_ba_2=sum(Ci_2.*Zi_2./Ai_2)/sum(Ci_2./Ai_2);
                    A_ba_2=1/sum(Ci_2./Ai_2);
                    J_ba_2=exp(sum(Ci_2.*Zi_2./Ai_2.*log(Ji_2))/sum(Ci_2.*Zi_2./Ai_2));
                    %����У�����dE/dS���Լ���������
                    del_2=-7.85*10^4*rho_2*Z_ba_2/A_ba_2/E*log(1.166*E/J_ba_2+0.0063*Z_ba_2+0.7334);
                    
                    %���ԭ��������֮����룬��lamdaӦ��һ��
                    l=sqrt(((x_last-x_now))^2+((y_last-y_now))^2+((z_last-z_now))^2);
                    

                    %����Ӵ�����ǰ�ĳ��ȣ����ݹ�ʽ�����һ�γ���
                    l1=(z_last-depth_add(max(count_2,count))*10^-7)/(z_last-z_now)*l;
                    lamda_2=(lamda-l1)*lamda_ba_2/lamda_ba;
                    %�����������ɢ�䲽��lamda
                    lamda= l1+lamda_2;
                    %�����������������ʧ��
                    del=(del*l1+del_2*lamda_2)/lamda;
                    
                    %����������ĵ���ɢ���յ������
                    x_now=(x_now-x_last)*lamda/l+x_last;
                    y_now=(y_now-y_last)*lamda/l+y_last;
                    z_now=(z_now-z_last)*lamda/l+z_last;
                        
                end
                break;
            end
        end  
        
        %���������������������del��lamda���������ڻ�ͼ
        E=E+del*lamda;
        del_energy=[del_energy,del];
        lamda_energy=[lamda_energy,lamda];
      
        %�����˱߽�����
        if z_now*10^7>0
            %������ڹ���֮�У����������
            finalx=[finalx,x_now*10^7];finaly=[finaly,y_now*10^7];finalz=[finalz,z_now*10^7];
        else
            %����Ѿ��ݳ�������棬���E��Ϊ�㣬��������Ϊ��ɢ�����
            E=0;
            is_back=1;
        end
    end
    %�����ӹ켣��ͼ
    figure(1);title('Electronic trajectory3D');xlabel('x/nm');ylabel('y/nm');zlabel('z/nm');
    %�ж��Ƿ�ɢ����ӣ��ò�ͬ��ɫ��ע
    if is_back==1
        plot3(finalx,finaly,-finalz,'r');hold on
    else
        plot3(finalx,finaly,-finalz,'b');hold on
    end   
    figure(5);title('Electronic trajectory');xlabel('x/nm');ylabel('z/nm');
    if is_back==1
        count_back=count_back+1;
        plot(finaly,-finalz,'r');hold on
    else
        plot(finaly,-finalz,'b ');hold on
    end 
    %���ش˵��ӵĺ�������ֵ
    miny=min(finaly);maxy=max(finaly);
    minx=min(finalx);maxx=max(finalx);

    %�����һ��ɢ���ĵ��������Լ�������ʧ�ʡ�������ʧ�Ա�������
    x_x=[x_x,x_now*10^7];
    y_x=[y_x,y_now*10^7];
    z_x=[z_x,z_now*10^7];
    del_energy=[del_energy,E];
    lamda_energy=[lamda_energy,-1];
    %���������ֲ�
    energy_res=energy_distribution(z_x,del_energy,lamda_energy,zmax,energy_res);
    %����������r��ʱ�������ͼ
    r_x=sqrt(x_x.^2+y_x.^2);
    if is_back==1
        %��ɢ�����
        energy_res_back=energy_distributionR(r_x,del_energy,lamda_energy,energy_res_back,z_x,layer{1,5},hengzuobiao);
    else
        %ǰɢ�����
        energy_res_forward=energy_distributionR(r_x,del_energy,lamda_energy,energy_res_forward,z_x,layer{1,5},hengzuobiao);
    end
    %���˴ε���ɢ������ݾ����浽һ��cell�У���һ������Ϊ����ɢ���x���꣬�ڶ���Ϊy��������Ϊz�����ĸ�Ϊÿ��ɢ����ʧ������
    electron_temp={};
    electron_temp{1}=finalx;
    electron_temp{2}=finaly;
    electron_temp{3}=finalz;
    electron_temp{4}=del.*lamda;
    %������������ֵ�����Է����ͼ
    coor_temp=[min(finalx),max(finalx),min(finaly),max(finaly),min(finalz),max(finalz)];
end

