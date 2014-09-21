clear all;close all;clc;tic;
%***************************
%整个算法流程还是比较简单，和传统的MUSCI算法差不多，你只需要写TCM和SCM的算法流程即可
%而这个两个算法你可以写在一个文件里，不同函数调用即可。

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
        %C++代码从此处开始写，需要写的部分不多，但是你首先必须用清楚，你需要哪些
        %变量从外部传进来。比如：func(int a, char b)，你得知道你需要变量a和b。
        %在这里你可以先写一个简单的程序：
        %spectrum_len：空间谱的长度，一般我们
        %设置的好像是1024，你还要知道信号的长度(也就是采样点数，程序里是snapshots)
        %其他可能会用到的变量，你到时可以问我，写的时候
        %你可以参考specest_music_armadillo_impl.cc
        
        %这里norm代表范数，norm(*,2)代表2范数，A‘代表共轭转置，其他还有不懂的
        %符号或者函数，你可以问我，整个程序你不懂也可以问我。
        
        %如果你想要测试你所写的程序是否正确，你可以将上面的X=A*S+noise这里的X读
        %入文件保存为dat文件，然后按照read_from_file里的程序写测试程序即可。
        %% covariance matrix
        R=X*X'/snapshots;%传统的music算法，得到信号的协方差
        %这里是共轭转置A.t() or trans(A)
        
        %这一块是SCM算法
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
        %这一块是TCM算法
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
        %下面的程序是对三种不同方法得到的协方差(R)，协变异函数(Cscm,Ctcm)
        %进行特征值分解
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
        %然后在此次得到空间谱
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