clear;
clc;
close all;
%生成仿真数据,8基因，20个样本，25个时间节点
e=exp(1);
for i=1:8  %network node
    X(i,1)=rand()*(5-1)+1;%生成初始值
end
delta_t=1;
L=2100;       
N=round(L/delta_t);  
ts=zeros(N,1);
sigma0=0.06;  %noise
es=0.000000001;
s0(1)=-0.5;
s0(2)=-0.475;
s0(3)=-0.45;
s0(4)=-0.4;
s0(5)=-0.375;
s0(6)=-0.35;
s0(7)=-0.325;
s0(8)=-0.3;
s0(9)=-0.275;
s0(10)=-0.25;
s0(11)=-0.225;
s0(12)=-0.2;
s0(13)=-0.175;
s0(14)=-0.15;
s0(15)=-0.125;
s0(16)=-0.1;
s0(17)=-0.08;
s0(18)=-0.05;
s0(19)=-0.02;
s0(20)=-0.001;%临界点
s0(21)=0.01;
s0(22)=0.05;
s0(23)=0.08;
s0(24)=0.1;
s0(25)=0.15;

D= [-1 1 0 0 0 0 0 0 ;...
    -1 -1 0 0 0 0 0 0 ;...
    1 0 -1 0 0 0 0 0 ;...
    1 0 0 -1 0 0 0 0 ;...
    1 0 0 0 -1 0 0 0 ;...
    0 0 0 0 1 -1 0 1 ;...
    0 0 0 0 0 0 -1 1 ;...
    0 0 0 0 0 0 -1 -1];
TT=100;
stage=25;
weight_SD=zeros(1,stage);
sample_num=1;%案例样本数
CC=zeros(8,4,sample_num);  %案例样本矩阵 
CC1=zeros(8,4,20);   %参考样本矩阵      
%for i=1:20
%    ss(i)=s(28-i);
%end
ss=-0.5:5/1900:-0.45;%定义状态变量范围
%1. reference samples
for ll=1:20
    qq(ll)=0.96^(1/abs(ss(ll)));
    E=[-3/5*qq(ll) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
    J=D*E*inv(D);
    for i=1:N-1
        ts(i+1)=ts(i)+delta_t;
        eJ=e^(J*delta_t);
        for jj=1:8
            X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma0*normrnd(0,1)*delta_t;
        end
    end
    CC1(:,1,ll)=X(:,2000); %存储样本数据    
end
pprofile=reshape(CC1(:,1,:),8,20);%8个基因，20个参考样本，25个时间节点
ssize=size(pprofile);
%2. case samples
HD=zeros(8,25);
network={};
fid=fopen("adj_edges_all.txt");
[network,total]=signet(fid);
mi=zeros(total,total,25);
for T=1:TT
tempcase=zeros(stage,8,21);% reference+1 case
patients_num=21;%病例数量
for t=1:25
    q(t)=0.96^(1/abs(s0(t)));%计算状态变量对应的q值
    E=[-3/5*q(t) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
    J=D*E*inv(D);
    for k=1:sample_num%案例样本
        for i=1:N-1
            ts(i+1)=ts(i)+delta_t;%时间更新
            eJ=e^(J*delta_t);
            for jj=1:8
                X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma0*normrnd(0,1)*delta_t;
            end
        end
        CC(:,t,k)=X(:,2000);%存储样本数据
    end
    TC(:,t)=reshape(CC(:,t,:),8,sample_num);
    tempcase(t,:,1:patients_num-1)=pprofile(:,1:ssize(2));%tempcase(时间节点，基因，样本），每个时间节点的参考样本相同
    tempcase(t,:,patients_num)=TC(:,t);
%     tempcase=abs(tempcase);%取绝对值
end

%stages
tempcase=permute(tempcase,[2,3,1]);%调整数组维度顺序，变成tempcase(基因，样本，时间节点)
psize=size(tempcase);%获取数组大小
case_num=20*ones(1,25);%每个时间节点的参考样本数量
bnd_cdf=zeros(psize(1),psize(1),25);  %初始化边界累积分布函数数组

% tempcase=log(1+tempcase);%对数组取对数并加一处理
MI=zeros(total,stage);%初始化互信息矩阵
for t=1:stage
    for c=1:psize(1)
        mu(c)=mean(tempcase(c,1:case_num(t),t));%样本均值
        sigma(c)=std(tempcase(c,1:case_num(t),t)); %样本标准差
    end
    for na=1:total        
        edges_list=[];% 存储与节点相邻的边的列表
        center=network{na}{1};% 获取当前节点的中心节点编号
        e=0;% 初始化边的数量为0
        for n=2:length(network{na})  % 遍历与当前节点相邻的节点
            nei=network{na}{n};  % 获取相邻节点的编号
            e=e+1;  % 边数加1
            edges_list(e, :)=[str2num(center) str2num(nei)];  % 添加边到边列表
        end
        if e<1
            continue; % 如果没有边，继续下一个节点
        end
        c=str2num(center);
        casecenter_cdf=normcdf(tempcase(c,patients_num,t),mu(c),sigma(c))+es; %计算累积分布函数值        
        for i = 1:e%局部网络中的边
            casenei_cdf=normcdf(tempcase(edges_list(i, 2),patients_num,t),mu(edges_list(i, 2)),sigma(edges_list(i, 2)))+es; %计算累积分布函数值
            if (tempcase(edges_list(i, 1),patients_num,t)==0)||(tempcase(edges_list(i, 2),patients_num,t)==0)
                bnd_cdf=0;
            else
                rou=corrcoef(tempcase(edges_list(i, 1),1:case_num(t),t),tempcase(edges_list(i, 2),1:case_num(t),t));%计算相关系数
                bnd_pdf=@(u,v)1./(2.*pi.*sigma(edges_list(i, 1)).*sigma(edges_list(i, 2)).*sqrt(1-rou(1,2).*rou(1,2))).*exp(-1./(2.*(1-rou(1,2).^2)).*[(u-mu(edges_list(i, 1))).*(u-mu(edges_list(i, 1)))./(sigma(edges_list(i, 1)).*sigma(edges_list(i, 1)))-2.*rou(1,2).*(u-mu(edges_list(i, 1))).*(v-mu(edges_list(i, 2)))./(sigma(edges_list(i, 1)).*sigma(edges_list(i, 2)))+(v-mu(edges_list(i, 2))).*(v-mu(edges_list(i, 2)))./(sigma(edges_list(i, 2)).*sigma(edges_list(i, 2)))]);%计算边界概率密度函数
                bnd_cdf=integral2(bnd_pdf,0,tempcase(edges_list(i, 1),patients_num,t),0,tempcase(edges_list(i, 2),patients_num,t))+es;%计算边界累积分布函数
            end
            if bnd_cdf~=0
                w=bnd_cdf*log(bnd_cdf/(casecenter_cdf*casenei_cdf));%计算互信息
                mi(edges_list(i, 1),edges_list(i, 2),t)=mi(edges_list(i, 1),edges_list(i, 2),t)+abs(w);
                MI(c,t)=MI(c,t)+abs(w);%累加
            end
        end
        MI(c,t)=MI(c,t)/e;
        %算SD得到一个8*25矩阵，和局部网络熵8*25矩阵对应相乘，得到一个8*25矩阵，每一列取前5个平均值作为指标得到1*25矩阵做曲线
        sd_n=std(tempcase(c,1:case_num(t),t));
        sd_add=std(tempcase(c,:,t));
        SD(c,t)=abs(sd_n-sd_add);
        disp(['Time' num2str(t) 'Gene' num2str(na) 'is completed']);%显示计算完成信息
    end
end
HD=sqrt(MI.*SD);
for t=1:stage
    [HD_sorted,~]=sort(HD(:,t),'descend');
    weight_SD(t)=weight_SD(t)+mean(HD_sorted(1:3));
end
end
mi=mi/TT;
weight_SD=weight_SD/TT;

%曲线图
plot(s0(1:25),weight_SD,'r-*','LineWidth',2.5);%绘制加权差分熵随参数p变化的图像
%ylim([0 4.5])
xlabel('Parameter p');%设置X轴标签
ylabel('Score');%设置Y轴标签