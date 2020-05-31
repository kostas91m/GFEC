clc;
clear;
X1 = dlmread("1xData.dat");
Y1 = dlmread("1yData.dat");
Z1 = dlmread("1zData.dat");
X2 = dlmread("2xData.dat");
Y2 = dlmread("2yData.dat");
Z2 = dlmread("2zData.dat");
hold on;
contourf(X1,Y1,Z1,20)
contourf(X2,Y2,Z2,20)
colorbar
%caxis([0 14])
colormap("jet")
hold off;