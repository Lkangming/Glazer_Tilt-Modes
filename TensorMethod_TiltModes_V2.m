%% 利用Deformation Tensor(two-order tensor)获得特定tilt mode(基于Glazer notation)下的单胞原子构型：

% data：2024.08.28 V1--->V2
% author：K.M.Luo
% 参考文献：Beanland R. Structure of planar defects in tilted perovskites. Acta Crystallogr A. 2011 

%% 基本假设：1）氧原子"角共享"；2）氧原子位移不改变氧八面体的空间反演对称性(也即单个氧八面体上仅3个独立的氧原子，共9个自由度需确定，故只需要一个二阶张量-deformation tensor即可；）

% input：tilt modes(Glazer notations)+POSCAR+特定原子idx(其中"！！！待输入"字样为脚本当前的必要输入参数)
% output：特定tilt mode下超胞中的各原子坐标
% 目前仅考虑单胞双倍增效应的畸变

clear all

%% 1.POSCAR 数据读取：

% 打开POSCAR文件
% filename = 'POSCAR.vasp'; %待输入的文件名
filename = 'POSCAR.vasp';
fid = fopen(filename, 'r');
% 读取系统名称
system_name = fgetl(fid);
% 读取缩放因子
scale_factor = str2double(fgetl(fid));
% 读取晶格矢量
lattice_vectors = zeros(3, 3);
for i = 1:3
    line = fgetl(fid);
    lattice_vectors(i, :) = scale_factor * sscanf(line, '%f %f %f');
end
% 读取元素信息
elements_line = fgetl(fid);
if ~isempty(str2num(elements_line))  % 判断是否是数值
    elements = {};  % 如果是数值，说明没有元素信息
    atomic_counts = str2num(elements_line);
else
    elements = strsplit(elements_line);  % 解析元素种类
    atomic_counts = str2num(fgetl(fid));  % 读取元素的原子个数
end
% 读取原子坐标类型（Direct或Cartesian）
coord_type = strtrim(fgetl(fid));
% 读取原子坐标
num_atoms = sum(atomic_counts);
atom_coords = zeros(num_atoms, 3);
for i = 1:num_atoms
    line = fgetl(fid);
    atom_coords(i, :) = sscanf(line, '%f %f %f');
end
% 关闭文件
fclose(fid);

%% 2.基本参数输入
alpha1 = pi/20; %！！！待输入
beta1 = pi/20;
gama1 = 0;

%tilt modes：
 %以a-a-c+为例：即表明t(000)=(a,b,c) t(111)=(-a,-b,c)
 %以a0a0c-为例：即表面t(000)=(0,0,c) t(111)=(0,0,-c)
   %(000)位置上的tilt vector：
    a1 = sin(alpha1); 
    b1 = sin(beta1);
    c1 = sin(gama1);
   %(111)位置上的tilt vector：
    a2=-a1; b2=-b1; c2=c1;  %！！！待输入(该关系可由输入的Glazer符号直接确定，待后续优化)

%deformation tensor的对角元
a0 = cos(beta1)+cos(gama1)-2;
b0 = cos(alpha1)+cos(gama1)-2;
c0 = cos(alpha1)+cos(beta1)-2;

%特征长度因子：
l1 = 7.9381/4; l2=7.9381/4; l3=7.9381/4; %！！！待输入，如是特殊的正交基矢则可直接读入，但若是非正交基则需根据特定几何关系写入

%局部原点：
order_ori = [1,5,7,3,2,6,8,4];%！！！待输入(八面体中心原子在POSCAR中的idx)，但应可由晶体学相关操作直接读出POSCAR中特定位置对应的序号——>todo
ORI = zeros(3,8);
idx = zeros(1,8);

for j = 1:length(order_ori)
    idx = order_ori(1,j);
    for i = 1:3
        ORI(i,j) = atom_coords(idx,i);
    end
end

%% 3.tilt vectors和deformation tensors的确定

%判据量(特定位置上Q矩阵与q向量的正负判断因子，详见论文)：
one = ones(3,1);
V = zeros(1,8);
for i = 1:8
    V(:,i) = (-1)^(one.'*((ORI(:,i)-ORI(:,1))./0.5));
end

 %operator matrix Q and q：
 Q0 = zeros(3,3); q0 = zeros(3,1);
 Q0(1,1) = (a1+a2)/l1; Q0(2,2) = (b1+b2)/l2; Q0(3,3) = (c1+c2)/l3;
 q0(1,1) = a1/l1; q0(2,1) = b1/l2; q0(3,1) = c1/l3;

  %特定位点上的operator matrix：
Q_Matrix = cell(1,8);
q_vector = zeros(3,8);
for i = 1:8
    Q_Matrix(1,i) = {Q0.*(-V(1,i))};
    q_vector(:,i) = q0.*V(1,i);
end

  %特定位点上的tilt vector：
Tilt_vector = zeros(3,8);
for i = 1:8
    Tilt_vector(:,i) = cell2mat(Q_Matrix(1,i))*((ORI(:,i)-ORI(:,1))./0.5)+q_vector(:,i);
end


  %特定位点上的deformation tensor：
D = cell(1,8);
for i = 1:8
    D(1,i) = {[a0,Tilt_vector(3,i),-Tilt_vector(2,i);-Tilt_vector(3,i),b0,Tilt_vector(1,i);Tilt_vector(2,i),-Tilt_vector(1,i),c0]};
end

%% 4.按给定tilt modes操纵原子坐标以实现倾转构型：

%Oxygen atoms的原坐标：
cord0 = cell(3,8);
order_oxgen = [21,17,19,23,22,18,20,24;...
               27,31,29,25,28,32,30,26;...
               34,38,40,36,33,37,39,35]; %！！！待输入，但应该可通过一些晶体学的方法直接从POSCAR中读取特定位置的原子idx？待后续优化
atom_coords_trans = atom_coords.';

for j = 1:8
    for i = 1:3
        cord0(i,j) = {atom_coords_trans(:,order_oxgen(i,j))};
    end
end

 
 %倾转后Oxygen atoms的坐标：
 cord1 = cell(3,8);
 for j = 1:8
     k = cell2mat(D(1,j));
     for i = 1:3
%          k(i,i)=0; %保约束
         cord1(i,j)={cell2mat(cord0(i,j))+k(:,i)/4};
     end
 end

%% 5.输出具有特定Glazer tilt modes结构的POSCAR


 Cord_list = [cell2mat(cord1(1,2)) cell2mat(cord1(1,6)) cell2mat(cord1(1,3)) cell2mat(cord1(1,7)) cell2mat(cord1(1,1)) cell2mat(cord1(1,5)) ...
     cell2mat(cord1(1,4)) cell2mat(cord1(1,8)) cell2mat(cord1(2,4)) cell2mat(cord1(2,8)) cell2mat(cord1(2,1)) cell2mat(cord1(2,5)) ...
     cell2mat(cord1(2,3)) cell2mat(cord1(2,7)) cell2mat(cord1(2,2)) cell2mat(cord1(2,6)) cell2mat(cord1(3,5)) cell2mat(cord1(3,1)) ...
     cell2mat(cord1(3,8)) cell2mat(cord1(3,4)) cell2mat(cord1(3,6)) cell2mat(cord1(3,2)) cell2mat(cord1(3,7)) cell2mat(cord1(3,3))];
%  outcome = Cord_list.'
Cord_list_trans = Cord_list.';

atom_coords_tilt = zeros(40,3);
for i = 1:40
    if i<=16
        atom_coords_tilt(i,:) = atom_coords(i,:);
    else
        atom_coords_tilt(i,:) = Cord_list_trans(i-16,:);
    end
end

% 写入POSCAR_tilt文件
title = 'Ti1 Pb1 O3';
Newfile = fopen('POSCAR_tilt_a-a-c0.vasp','w');
 %写入相关数据：
 fprintf(Newfile, '%s\n',title); %标题
 fprintf(Newfile, '%f\n',scale_factor); %缩放因子
 for i = 1:size(lattice_vectors,1) %正空间格矢
     fprintf(Newfile, '%f %f %f\n',lattice_vectors(i,:));
 end
 fprintf(Newfile, '%s ', elements{:});%元素符号
 fprintf(Newfile,'\n'); 
 fprintf(Newfile, '%d ',atomic_counts);%元素数目
 fprintf(Newfile,'\n');
 fprintf(Newfile, '%s\n',coord_type);%坐标类型
 for i = 1:size(atom_coords_tilt,1)%原子坐标录入
     fprintf(Newfile, '%f %f %f\n',atom_coords_tilt(i,:));
 end

 fclose(Newfile);
















