function [a,t]=signet(fid)
%���ļ�fid��ÿһ�ж�����һ�����磬������a;t��ʾ�������
adjacent_network={};
j=0;
while ~feof(fid)
    tline=fgetl(fid);
    j=j+1;
    adjacent_network{j}=regexp(tline, '\t', 'split');
end
fclose(fid);
t=j;a=adjacent_network;
