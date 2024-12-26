clear;
clc;
close all;
%%
%导入数据
data=xlsread("TCGA_LIHC.xlsx");
normal=data(:,1:50);%正常样本
case_mprofile=data(:,51:360);%疾病样本

network={};
fid=fopen("Gene_network.txt");
[network,total]=signet(fid);

reference_sample_num = 50;% 基准样本数量
patients_num=[157,78,70,5];% 各阶段病例数量
es=0.000000001;

tempcase(:, 1, 1:patients_num(1)) = case_mprofile(:, 1:157);   % Stage I
tempcase(:, 2, 1:patients_num(2)) = case_mprofile(:, 158:235); % Stage II
tempcase(:, 3, 1:patients_num(3)) = case_mprofile(:, 236:305); % Stage III
tempcase(:, 4, 1:patients_num(4)) = case_mprofile(:, 306:310); % Stage IV


psize=size(tempcase);%获取数组大小
HD= zeros(total,157,4);
weight_SD=zeros(1,4);
%tempcase=log(1+tempcase);%对数组取对数并加一处理

%%

mu=mean(normal(:,:)');%样本均值
sigma=std(normal(:,:)'); %样本标准差
for l=1:psize(2)% 遍历各个时间节点
    for s=1:patients_num(l)%遍历各个时间节点的病例样本
        for na=1:total%遍历每个网络节点
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

            clear bnd_cdf
            casecenter_cdf=normcdf(tempcase(str2num(center),l,s),mu(str2num(center)),sigma(str2num(center)))+es; %计算累积分布函数值
            MI=0;
            for i = 1:e%局部网络中的边
                casenei_cdf=normcdf(tempcase(edges_list(i, 2),l,s),mu(edges_list(i, 2)),sigma(edges_list(i, 2)))+es; %计算累积分布函数值
                if (tempcase(edges_list(i, 1),l,s)==0)||(tempcase(edges_list(i, 2),l,s)==0)
                    bnd_cdf=0;
                else
                    rou=corrcoef(normal(edges_list(i, 1),:),normal(edges_list(i, 2),:));%计算相关系数
                    bnd_pdf=@(u,v)1./(2.*pi.*sigma(edges_list(i, 1)).*sigma(edges_list(i, 2)).*sqrt(1-rou(1,2).*rou(1,2))).*exp(-1./(2.*(1-rou(1,2).^2)).*[(u-mu(edges_list(i, 1))).*(u-mu(edges_list(i, 1)))./(sigma(edges_list(i, 1)).*sigma(edges_list(i, 1)))-2.*rou(1,2).*(u-mu(edges_list(i, 1))).*(v-mu(edges_list(i, 2)))./(sigma(edges_list(i, 1)).*sigma(edges_list(i, 2)))+(v-mu(edges_list(i, 2))).*(v-mu(edges_list(i, 2)))./(sigma(edges_list(i, 2)).*sigma(edges_list(i, 2)))]);%计算边界概率密度函数
                    bnd_cdf=integral2(bnd_pdf,0,tempcase(edges_list(i, 1),l,s),0,tempcase(edges_list(i, 2),l,s))+es;%计算边界累积分布函数
                end
                if bnd_cdf~=0
                    mi=bnd_cdf.*log(bnd_cdf./(casecenter_cdf.*casenei_cdf));%计算互信息                    
                    MI=MI+abs(mi);%累加
                end
            end
            MI=MI/e;
            sd_n=std(normal(str2num(center),:));
            sd_add=std([normal(str2num(center),:),reshape(tempcase(str2num(center),l,s),1,1)]);
            SD=abs(sd_n-sd_add);
            HD(na,s,l)=MI*SD;
    end
    disp(['Time' num2str(l) 'Sample' num2str(s) ' is completed']);%显示计算完成信息
    end
end
HD_size = size(HD);
case_result = zeros(157, 4);
% 计算每个阶段的HD的平均值
for t = 1:HD_size(3)
    for case_num = 1:patients_num(t)       
        [sort_HD, idx] = sort(HD(:, case_num, t), 'descend');
        case_result(case_num, t) = mean(sort_HD(1:400));
    end
    result(t) = mean(case_result(1:patients_num(t), t));
end
%%
% 画图
plot(1:4,result,'r-*','LineWidth',2.5);%绘制加权差分熵随参数p变化的图像
B={'I' 'II' 'III' 'IV'};
set(gca,'XTick',1:4);
set(gca,'XTickLabel',B, 'FontWeight', 'bold', 'FontSize', 14);
xlabel('Stage', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('Score', 'FontWeight', 'bold', 'FontSize', 14);
line([2 2],[0.22 result(2)],'linestyle','--','Color','m','LineWidth',2);
%% 
% 筛选信号基因
%%%%%%%%% selected marker genes %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 定义熵的矩阵大小
HD_size = size(HD);  

% 初始化变量
selected_deltaH_genes = [];
fid = fopen('selected_LIHC_genes.txt', 'w'); 
selected_deltaH_genes = containers.Map;

% 对第二类病人进行操作
for s = 1:patient_num(2)  
    [tmp_com_idx, idx] = sort(HD(:, s, 2), 'descend');  

    % 选择熵值最高的5%的基因
    tmp_genes = node1(idx(1:floor(HD_size(1) * 0.05))); 

    % 计数每个选中的基因出现的次数
    for i = 1:length(tmp_genes)  % 遍历选中的基因
        if isKey(selected_deltaH_genes, tmp_genes{i}) == 0 
            selected_deltaH_genes(tmp_genes{i}) = 1;
        else
            selected_deltaH_genes(tmp_genes{i}) = selected_deltaH_genes(tmp_genes{i}) + 1;  
        end
    end
end

% 获取所有选中的基因，并写入文件
genes = selected_deltaH_genes.keys();  
for i = 1:length(genes) 
    fprintf(fid, '%s\t%d\n', genes{i}, selected_deltaH_genes(genes{i})); 
end
fclose(fid);  % 关闭文件
