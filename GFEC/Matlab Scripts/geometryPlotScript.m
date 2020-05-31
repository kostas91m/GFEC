clc;
clear;
nodeCoord = dlmread("coordinateData.dat");
connect = dlmread("connectivityData.dat");
contactConnect = dlmread("contactConnectivityData.dat");

for i=1:1:112
  node1 = connect(i,2);
  node2 = connect(i,3);
  node3 = connect(i,4);
  node4 = connect(i,5);
  node1XY = [nodeCoord(node1, 2), nodeCoord(node1, 3)];
  node2XY = [nodeCoord(node2, 2), nodeCoord(node2, 3)];
  node3XY = [nodeCoord(node3, 2), nodeCoord(node3, 3)];
  node4XY = [nodeCoord(node4, 2), nodeCoord(node4, 3)];
  
  nodesX = [nodeCoord(node1, 2) nodeCoord(node2, 2) nodeCoord(node3, 2) nodeCoord(node4, 2) nodeCoord(node1, 2)];
  nodesY = [nodeCoord(node1, 3) nodeCoord(node2, 3) nodeCoord(node3, 3) nodeCoord(node4, 3) nodeCoord(node1, 3)];
  hold on;
  line(nodesX, nodesY, 'LineWidth', 2, 'Color', 'blue');
  %plot(node2XY, node3XY);
  %plot(node3XY, node4XY);
  %plot(node4XY, node1XY);
endfor

for i=1:1:8
  node1 = contactConnect(i,2);
  node2 = contactConnect(i,3);
  
  node1XY = [nodeCoord(node1, 2), nodeCoord(node1, 3)];
  node2XY = [nodeCoord(node2, 2), nodeCoord(node2, 3)];
 
  
  nodesX = [nodeCoord(node1, 2) nodeCoord(node2, 2) ];
  nodesY = [nodeCoord(node1, 3) nodeCoord(node2, 3) ];
  hold on;
  p2 = line(nodesX, nodesY, 'LineWidth', 2, 'Color', 'red');
  %plot(node2XY, node3XY);
  %plot(node3XY, node4XY);
  %plot(node4XY, node1XY);
endfor
