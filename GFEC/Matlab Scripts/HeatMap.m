clc;
clear;
totalNodes = 486;
nodesInXcoor = 81;
nodesInYcoor = 3;
data = dlmread("Results30.dat");
data(1,:)=[]

upperBeamData = data(1:totalNodes/2, 1);

xDataUpper = data(1:totalNodes/2,1);
yDataUpper = data(1:totalNodes/2,2);
tempDataUpper = data(1:totalNodes/2,3);

xDataLower = data(totalNodes/2+1:totalNodes,1);
yDataLower = data(totalNodes/2+1:totalNodes,2);
tempDataLower = data(totalNodes/2+1:totalNodes,3);


contourXDataUpper=reshape(xDataUpper,nodesInXcoor,nodesInYcoor);
contourYDataUpper=reshape(yDataUpper,nodesInXcoor,nodesInYcoor);
contourtempDataUpper=reshape(tempDataUpper,nodesInXcoor,nodesInYcoor);

contourXDataLower=reshape(xDataLower,nodesInXcoor,nodesInYcoor);
contourYDataLower=reshape(yDataLower,nodesInXcoor,nodesInYcoor);
contourtempDataLower=reshape(tempDataLower,nodesInXcoor,nodesInYcoor);


%contour(xData, yData);
hold on;
contourf(contourXDataUpper, contourYDataUpper, contourtempDataUpper);
contourf(contourXDataLower, contourYDataLower, contourtempDataLower);
