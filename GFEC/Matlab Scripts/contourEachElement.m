clc;
clear;
connect = dlmread("connectivityData.dat");
finalSolutionData = dlmread("Results5.dat");

for i=1:1:320
  connect
  node1 = connect(i,2);
  node2 = connect(i,3);
  node3 = connect(i,4);
  node4 = connect(i,5);
  
  tableX = [finalSolutionData(node4, 1), finalSolutionData(node3, 1) ;
            finalSolutionData(node1, 1), finalSolutionData(node2, 1)];
            
  tableY = [finalSolutionData(node4, 2), finalSolutionData(node3, 2) ;
            finalSolutionData(node1, 2), finalSolutionData(node2, 2)];
            
  tableTemp = [finalSolutionData(node4, 3), finalSolutionData(node3, 3) ;
               finalSolutionData(node1, 3), finalSolutionData(node2, 3)];
  
  contourf(tableX, tableY, tableTemp);         
  hold on;
endfor


