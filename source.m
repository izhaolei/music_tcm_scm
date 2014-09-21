clear all;close all;clc;tic;
%***************************
%�����㷨���̻��ǱȽϼ򵥣��ʹ�ͳ��MUSCI�㷨��࣬��ֻ��ҪдTCM��SCM���㷨���̼���
%����������㷨�����д��һ���ļ����ͬ�������ü��ɡ�

%% paramter
J=sqrt(-1);sensor_num=4;lambda=1;d=lambda/2;snapshots=300;
src_num=2;theta=[0,10]*pi/180;
gamma=1;alpha_set=[2,1.4,1];beta=0;delta=0;GSNR=10;

trials=2;
est_theta=-90:90;
spectrum_len=length(est_theta);
aver_R_spectrum=zeros(length(alpha_set),spectrum_len);
aver_SCM_spectrum=zeros(length(alpha_set),spectrum_len);
aver_TCM_spectrum=zeros(length(alpha_set),spectrum_len);
for ii=1:length(alpha_set)
    alpha=alpha_set(ii);
    
    R_spectrum=zeros(trials,spectrum_len);
    SCM_spectrum=zeros(trials,spectrum_len);
    TCM_spectrum=zeros(trials,spectrum_len);
    for jj=1:trials
        %% alpha stable noise
        noise=(stblrnd(alpha,beta,gamma,delta,sensor_num,snapshots)+...
            J*stblrnd(alpha,beta,gamma,delta,sensor_num,snapshots))/sqrt(2);
        p_n=sum(diag((noise*noise')/snapshots))/sensor_num;
        P_n=sqrt(p_n);
        Ps=(P_n^alpha)*10^(GSNR/10);
        %Ps=(GSNR/10)*log(10)/log(P_n);
        %% signal
        fs=200;
        t=0:1/fs:(snapshots-1)/fs;
        S=zeros(src_num,snapshots);
        for i=1:src_num
            f=rand()*100;
            S(i,:)=sqrt(Ps)*exp(J*2*pi*f*t);
        end
        
        %% steering vector
        A=zeros(sensor_num,src_num);
        for i=1:src_num
            A(:,i)=exp(-J*2*pi*d*(0:sensor_num-1)'*sin(theta(i))/lambda);
        end
        
        %% receive signal
        X=A*S+noise;
        %C++����Ӵ˴���ʼд����Ҫд�Ĳ��ֲ��࣬���������ȱ��������������Ҫ��Щ
        %�������ⲿ�����������磺func(int a, char b)�����֪������Ҫ����a��b��
        %�������������дһ���򵥵ĳ���
        %spectrum_len���ռ��׵ĳ��ȣ�һ������
        %���õĺ�����1024���㻹Ҫ֪���źŵĳ���(Ҳ���ǲ�����������������snapshots)
        %�������ܻ��õ��ı������㵽ʱ�������ң�д��ʱ��
        %����Բο�specest_music_armadillo_impl.cc
        
        %����norm��������norm(*,2)����2������A��������ת�ã��������в�����
        %���Ż��ߺ�������������ң����������㲻��Ҳ�������ҡ�
        
        %�������Ҫ��������д�ĳ����Ƿ���ȷ������Խ������X=A*S+noise�����X��
        %���ļ�����Ϊdat�ļ���Ȼ����read_from_file��ĳ���д���Գ��򼴿ɡ�
        %% covariance matrix
        R=X*X'/snapshots;%��ͳ��music�㷨���õ��źŵ�Э����
        %�����ǹ���ת��A.t() or trans(A)
        
        %��һ����SCM�㷨
        %%  spatial sign covariance matrix(SCM)
        Total_Sscm=0;
        for i=1:snapshots
            if X(:,i)==0%	 A.col(k)
                Total_Sscm=Total_Sscm+0;
            else
                Sscm=X(:,i)/norm(X(:,i),2);%norm( X, p )
                Total_Sscm=Total_Sscm+Sscm*Sscm';
            end
        end
        Cscm=Total_Sscm/snapshots;
        %��һ����TCM�㷨
        %% spatial Kendall's tau covariance matrix(TCM)
        Total_Stcm=0;
        for i=1:snapshots
            for j=1:snapshots
                Xtcm=X(:,i)-X(:,j);
                if Xtcm==0
                    Total_Stcm=Total_Stcm+0;
                    %             fprintf('Xtcm==0\n')
                else
                    Stcm=Xtcm/norm(Xtcm,2);
                    Total_Stcm=Total_Stcm+Stcm*Stcm';
                    %             fprintf('Xtcm!=0\n')
                end
            end
        end
        Ctcm=Total_Stcm/(snapshots*(snapshots-1));
        %����ĳ����Ƕ����ֲ�ͬ�����õ���Э����(R)��Э���캯��(Cscm,Ctcm)
        %��������ֵ�ֽ�
        %% eig
        [R_V, R_D]=eig(R);
        [~,R_Ind]=sort(diag(R_D),'descend');
        R_Gn=R_V(:,R_Ind(src_num+1:end));
        
        [SCM_V, SCM_D]=eig(Cscm);
        [~,SCM_Ind]=sort(diag(SCM_D),'descend');
        SCM_Gn=SCM_V(:,SCM_Ind(src_num+1:end));
        
        [TCM_V, TCM_D]=eig(Ctcm);
        [~,TCM_Ind]=sort(diag(SCM_D),'descend');
        TCM_Gn=TCM_V(:,TCM_Ind(src_num+1:end));
        %Ȼ���ڴ˴εõ��ռ���
        %% music
        for i=1:spectrum_len
            a_theta=exp(-J*2*pi*d*(0:sensor_num-1)'*sin(est_theta(i)*pi/180)/lambda);
            R_spectrum(jj,i)=1/abs(a_theta'*R_Gn*R_Gn'*a_theta);
        end
        
        for i=1:spectrum_len
            SCM_a_theta=exp(-J*2*pi*d*(0:sensor_num-1)'*sin(est_theta(i)*pi/180)/lambda);
            SCM_spectrum(jj,i)=1/abs(SCM_a_theta'*SCM_Gn*SCM_Gn'*SCM_a_theta);
        end
        
        for i=1:spectrum_len
            TCM_a_theta=exp(-J*2*pi*d*(0:sensor_num-1)'*sin(est_theta(i)*pi/180)/lambda);
            TCM_spectrum(jj,i)=1/abs(TCM_a_theta'*TCM_Gn*TCM_Gn'*TCM_a_theta);
        end
        
    end
    
    aver_R_spectrum(ii,:)=sum(R_spectrum)/trials;
    aver_SCM_spectrum(ii,:)=sum(SCM_spectrum)/trials;
    aver_TCM_spectrum(ii,:)=sum(TCM_spectrum)/trials;
    
end

%% figure
figure()
subplot(331)
plot(est_theta,10*log10(aver_R_spectrum(1,:)),'r','linewidth',1)
hold on
line([0 0],[-20 40],'Color','k','linestyle',':')
hold on
line([10 10],[-20 40],'Color','k','linestyle',':')
hold off
xlabel('angle')
ylabel('amplitude(dB)')
title('MUSIC,alpha=2')

subplot(332)
plot(est_theta,10*log10(aver_TCM_spectrum(1,:)),'r','linewidth',1)
hold on
line([0 0],[-20 40],'Color','k','linestyle',':')
hold on
line([10 10],[-20 40],'Color','k','linestyle',':')
hold off
xlabel('angle')
ylabel('amplitude(dB)')
title('TCM-MUSIC,alpha=2')


subplot(333)
plot(est_theta,10*log10(aver_SCM_spectrum(1,:)),'r','linewidth',1)
hold on
line([0 0],[-20 40],'Color','k','linestyle',':')
hold on
line([10 10],[-20 40],'Color','k','linestyle',':')
hold off
xlabel('angle')
ylabel('amplitude(dB)')
title('SCM-MUSIC,alpha=2')

subplot(334)
plot(est_theta,10*log10(aver_R_spectrum(2,:)),'r','linewidth',1)
hold on
line([0 0],[-20 40],'Color','k','linestyle',':')
hold on
line([10 10],[-20 40],'Color','k','linestyle',':')
hold off
xlabel('angle')
ylabel('amplitude(dB)')
title('MUSIC,alpha=1.4')

subplot(335)
plot(est_theta,10*log10(aver_TCM_spectrum(2,:)),'r','linewidth',1)
hold on
line([0 0],[-20 40],'Color','k','linestyle',':')
hold on
line([10 10],[-20 40],'Color','k','linestyle',':')
hold off
xlabel('angle')
ylabel('amplitude(dB)')
title('TCM-MUSIC,alpha=1.4')

subplot(336)
plot(est_theta,10*log10(aver_SCM_spectrum(2,:)),'r','linewidth',1)
hold on
line([0 0],[-20 40],'Color','k','linestyle',':')
hold on
line([10 10],[-20 40],'Color','k','linestyle',':')
hold off
xlabel('angle')
ylabel('amplitude(dB)')
title('SCM-MUSIC,alpha=1.4')

subplot(337)
plot(est_theta,10*log10(aver_R_spectrum(3,:)),'r','linewidth',1)
hold on
line([0 0],[-20 40],'Color','k','linestyle',':')
hold on
line([10 10],[-20 40],'Color','k','linestyle',':')
hold off
xlabel('angle')
ylabel('amplitude(dB)')
title('MUSIC,alpha=1')

subplot(338)
plot(est_theta,10*log10(aver_TCM_spectrum(3,:)),'r','linewidth',1)
hold on
line([0 0],[-20 40],'Color','k','linestyle',':')
hold on
line([10 10],[-20 40],'Color','k','linestyle',':')
hold off
xlabel('angle')
ylabel('amplitude(dB)')
title('TCM-MUSIC,alpha=1')

subplot(339)
plot(est_theta,10*log10(aver_SCM_spectrum(3,:)),'r','linewidth',1)
hold on
line([0 0],[-20 40],'Color','k','linestyle',':')
hold on
line([10 10],[-20 40],'Color','k','linestyle',':')
hold off
xlabel('angle')
ylabel('amplitude(dB)')
title('SCM-MUSIC,alpha=1')
axis tight
toc;