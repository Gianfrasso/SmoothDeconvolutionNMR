% FvdB 160605
clear all; clf; clc;
set(gcf,'PaperPositionMode','auto');

load NMR_POTATO_DM.mat
X = DataSet.data(:,1:end-2);
Y = DataSet.data(:,end-1:end);
ObjLab = DataSet.Labels{1};
Time = char(DataSet.Labels{2});
Time = str2num(Time(1:end-2,:))';

index = 1:2^11;
X = X(:,index);
Time = Time(index);

clear DataSet

save go01