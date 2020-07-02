clc;
clear;

data1 = dlmread("ContactForces19.dat");
data2= dlmread("ContactForces20.dat");
data3 = dlmread("contactivity19.dat");
data4 = dlmread("contactivity20.dat");


contactivity1 = sum(data3);
contactivity2 = sum(data4);

for i=1:1:9
  data1Norm(i,1) = norm(data1(:,i));
endfor
