clc;
clear;
dataFiles = 20;
contactElements = 9;

%Read contact forces and contactivity .dat files
%-------------------------------------------------------------------------------
for i=1:1:dataFiles
  path = strcat("ContactForces",num2str(i),".dat");
  contactForceData{i,1} = dlmread(path);
endfor

for i=1:1:dataFiles
  path = strcat("contactivity",num2str(i),".dat");
  contactivityData{i,1} = dlmread(path);
endfor

%Get contact forces norm for each element per load step. Get total contact force
%per load step
%-------------------------------------------------------------------------------
for j=1:1:dataFiles
  for i=1:1:contactElements
    contactForcesNorms{j,1}(i,1) = norm(contactForceData{j,1}(:,i));
  endfor
endfor

for i=1:1:dataFiles
  totalContactForce(i,1) = sum(contactForcesNorms{i,1});  
endfor

%Get total contactivity per load step
%-------------------------------------------------------------------------------
for i=1:1:dataFiles
  contactivity(i,1) = sum(contactivityData{i,1});
endfor

%Plot totalContactForce with contactivity
%-------------------------------------------------------------------------------
plot(totalContactForce, contactivity);
