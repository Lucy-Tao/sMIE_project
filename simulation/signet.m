function [a,t]=signet(fid)
%将文件fid中每一行都当做一个网络，保存至a;t表示网络个数
adjacent_network={};
j=0;
while ~feof(fid)
    tline=fgetl(fid);
    j=j+1;
    adjacent_network{j}=regexp(tline, '\t', 'split');
end
fclose(fid);
t=j;a=adjacent_network;
