function patran = pat_pottery()
%
%  Manipulate ship hull shape
%
% Arthur: Yijun
%       07/06/2016

%% Initialization
clear vars
addpath(cd);
[FileName,PathName,FilterIndex] = uigetfile('.pat');
file = [PathName FileName];
p = patran_read_pat(file);

%% Input
cd(PathName);
Name = 'LNG_L_Full.pat';

% Step 0: scale main dimension
Old_dim = [274, 44.2 ,11]; % [Lpp, B, T]
New_dim = [266,41.6,11.2];
sca_ratio = New_dim./Old_dim;

% Step 1: extend parallel body
factor_ext = 1.7; % TO DO!!!
Xrange = [-42, 42];

% Step 2: broaden center flat body [* optional]
% Step 3: Change depth [* done in Torben's script]

p2 = p;
for i = 1:3
    p2.crd(:,i) = p.crd(:,i).*sca_ratio(i);
end
%% Visual
fig = figure();
x = p2.crd(:,1);
y = p2.crd(:,2);
z = p2.crd(:,3);

plot3(x,y,z,'r.')
axis equal

%% Parallel body - x dir
x2 = [];
y2 = [];
z2 = [];

len_sec_old = [Xrange(1) - min(x), Xrange(2)-Xrange(1), max(x) - Xrange(2)];
len_sec_new = [(Xrange(1) - min(x)) - 0.5*(factor_ext - 1)*(Xrange(2)-Xrange(1)),...
    (Xrange(2)-Xrange(1))*factor_ext, max(x) - Xrange(2) - 0.5*(factor_ext - ...
    1)*(Xrange(2)-Xrange(1))];
sca_ratio_L = len_sec_new./len_sec_old;

index_aft = find(x<Xrange(1));
index_bow = find(x>Xrange(2));

index_tmp_par1 = find(x >= Xrange(1));
index_tmp_par2 = find(x <= Xrange(2));
index_par = intersect(index_tmp_par1,index_tmp_par2);

% size(index_aft)
% size(index_par)
% size(index_bow)

flag_aft(1) = min(x(index_aft));
flag_aft(2) = max(x(index_aft));

flag_bow(1) = min(x(index_bow));
flag_bow(2) = max(x(index_bow));

flag_par(1) = min(x(index_par));
flag_par(2) = max(x(index_par));

% % aft
tmp_aft = x(index_aft);
for i = 1: length(tmp_aft)
    x2(index_aft(i)) = (tmp_aft(i)- flag_aft(1)).*sca_ratio_L(1)+flag_aft(1);
end

% % test
% tmp1 = x2;
% tmp1(tmp1==0) = -70;
% max(tmp1)
% min(tmp1)
% max(tmp1) - min(tmp1)

% % parallel body
x2(index_par) = x(index_par).*sca_ratio_L(2);

% % box
tmp_bow = x(index_bow);
for i = 1: length(tmp_bow)
    x2(index_bow(i)) = flag_bow(2) - (-tmp_bow(i)+ flag_bow(2)).*sca_ratio_L(3);
end

% % test
% tmp1 = x2;
% tmp1(tmp1==0) = 90;
% max(tmp1)
% min(tmp1)
% max(tmp1) - min(tmp1)

%% Central plane expansion - y dir


%% More freeboard - z dir



%% Visual
hold on
plot3(x2,y,z,'g+')

%% Write into patran
patran = p2;
patran.crd(:,1) = x2;
patran.crd(:,2) = y;
patran.crd(:,3) = z;

patran_write(patran,PathName,Name);