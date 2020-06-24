using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Xml.Schema;
using OpenTK.Input;
using System.Windows.Media.Animation;

namespace GFEC
{
    public static class ExportToFile
    {
        public static void ExportExplicitResults(Dictionary<int, double[]> solution, Dictionary<int, double> timeAtEachStep, int dofNumber, int intervals)
        {
            string[] lines = new string[solution.Count / intervals];
            int step = 0;
            int line = 0;
            while (step < solution.Count - 1)
            {
                double[] sol = solution[step];
                lines[line] = line.ToString() + " " + timeAtEachStep[step].ToString() + " " + sol[dofNumber].ToString();
                line = line + 1;
                step = step + intervals;
                //if (step >= solution.Count-1)
                //{
                //    break;
                //}
            }
            //File.WriteAllLines(@"D:\WriteLines2.txt", lines);
        }

        public static void ExportGeometryDataWithTemperatures(Dictionary<int, INode> nodesList, double[] temperatures, string path)
        {
            if (nodesList.Count != temperatures.Length)
            {
                throw new Exception("Mismatch beetween total nodes and temperatures vector");
            }
            
            string[] lines = new string[nodesList.Count+1];
            lines[0] = "X\tY\tTemperature";
            for (int i = 1; i <= nodesList.Count; i++)
            {
                lines[i] = nodesList[i].XCoordinate.ToString() + "\t" + nodesList[i].YCoordinate.ToString() + "\t" + temperatures[i - 1];
            }
            File.WriteAllLines(path, lines);
        }

        public static void ExportGeometryDataWithTemperatures(ISolver structuralSolution, Dictionary<int, double[]> thermalSolutionsList, int[] thermalBoundaryConditions)
        {
            Dictionary<int, INode> initialNodesList = structuralSolution.AssemblyData.Nodes;
            Dictionary<int, double[]> temperaturesForAllSteps = BoundaryConditionsImposition.CreateFullVectorListFromReduced(thermalSolutionsList, thermalBoundaryConditions);
            Dictionary<int, double[]> displacementsForAllSteps = structuralSolution.GetAllStepsFullSolutionVectors();
            foreach (KeyValuePair<int, double[]> vector in temperaturesForAllSteps)
            {
                Dictionary<int, INode> finalNodesList = Assembly.CalculateFinalNodalCoordinates(initialNodesList, displacementsForAllSteps[vector.Key]);
                double[] temperatures = vector.Value;
                string[] lines = new string[finalNodesList.Count + 1];
                lines[0] = "X\tY\tTemperature";
                if (finalNodesList.Count != temperatures.Length)
                {
                    throw new Exception("Mismatch beetween total nodes and temperatures vector");
                }
                for (int j = 1; j <= finalNodesList.Count; j++)
                {
                    lines[j] = finalNodesList[j].XCoordinate.ToString() + "\t" + finalNodesList[j].YCoordinate.ToString() + "\t" + temperatures[j - 1];
                }
                File.WriteAllLines("Results/Results"+vector.Key.ToString() + ".dat", lines);
            }  
        }

        public static void CreateContourDataForMatlab(double[] x, double[] y, double[] z, int rows, int columns, string path)
        {
            double[,] xContour = new double[rows, columns];
            double[,] yContour = new double[rows, columns];
            double[,] zContour = new double[rows, columns];

            //int[] arr = Enumerable.Repeat(42, 10000).ToArray();

            //string[] yData = new string[rows];
            //string temp;
            //for (int i = 0; i < rows; i++)
            //{
            //    for (int j = 0; j < columns; j++)
            //    {
            //        yContour[i, j] = y[i *columns+ j]; 
            //    }
            //}
            xContour = VectorOperations.ConvertVectorToMatrix(x, rows, columns);
            yContour = VectorOperations.ConvertVectorToMatrix(y, rows, columns);
            zContour = VectorOperations.ConvertVectorToMatrix(z, rows, columns);

            //for (int i = 0; i < rows; i++)
            //{
            //    for (int j = 0; j < columns; j++)
            //    {
            //        yData[i] = yData[i] + "\t" + yContour[i, j];
            //    }
            //}
            //File.WriteAllLines(@"C:\Users\Public\Documents\ContourDataY.dat", yData);

            MatrixOperations.PrintMatrixToFile(xContour, path + "xData.dat");
            MatrixOperations.PrintMatrixToFile(yContour, path + "yData.dat");
            MatrixOperations.PrintMatrixToFile(zContour, path + "zData.dat");

        }

        public static void ExportMatlabInitialGeometry(IAssembly assembly)
        {
            List<string> coordinateData = new List<string>();
            List<string> connectivityData = new List<string>();
            List<string> contactConnectivityData = new List<string>();

            foreach (var node in assembly.Nodes)
            {
                coordinateData.Add(node.Key.ToString() + "\t" + node.Value.XCoordinate.ToString() + "\t" + node.Value.YCoordinate.ToString());
            }

            foreach (var element in assembly.ElementsAssembly)
            {
                string line = element.Key.ToString();

                if (element.Value is ContactNtN2D)
                {
                    foreach (var elementNode in assembly.ElementsConnectivity[element.Key])
                    {
                        line = line + "\t" + elementNode.Value.ToString();
                    }
                    contactConnectivityData.Add(line);
                }
                else
                {
                    foreach (var elementNode in assembly.ElementsConnectivity[element.Key])
                    {
                        line = line + "\t" + elementNode.Value.ToString();
                    }
                    connectivityData.Add(line);
                }
            }

            File.WriteAllLines(@"C:\Users\Public\Documents\coordinateData.dat", coordinateData);
            File.WriteAllLines(@"C:\Users\Public\Documents\connectivityData.dat", connectivityData);
            File.WriteAllLines(@"C:\Users\Public\Documents\contactConnectivityData.dat", contactConnectivityData);
        }

        public static void ExportMatlabFinalGeometry(IAssembly assembly, double[] displacements)
        {
            List<string> finalCoordinateData = new List<string>();

            Dictionary<int, INode> nodesAtFinalState = Assembly.CalculateFinalNodalCoordinates(assembly.Nodes, displacements);

            foreach (var node in nodesAtFinalState)
            {
                finalCoordinateData.Add(node.Key.ToString() + "\t" + node.Value.XCoordinate.ToString() + "\t" + node.Value.YCoordinate.ToString());
            }
            File.WriteAllLines(@"C:\Users\Public\Documents\finalCoordinateData.dat", finalCoordinateData);
        }

        public static void ExportCondactivityForAllLoadSteps(List<Dictionary<int, double>> contactContactivityForEachStep)
        {
            int k = 0;
            foreach (var loadStep in contactContactivityForEachStep)
            {
                k = k + 1;
                double[] contactContactivity = loadStep.Values.ToArray();
                VectorOperations.PrintVectorToFile(contactContactivity, "Results/contactivity" + k.ToString() + ".dat");
            }
        }

        public static void ExportContactForcesForAllLoadSteps(Dictionary<int, Dictionary<int, double[]>> allStepsContactForces)
        {
            int k = 0;
            foreach (KeyValuePair<int, Dictionary<int, double[]>> loadStep in allStepsContactForces)
            {
                k = k + 1;
                double[] contactContactivity = loadStep.Values.ToArray();
                VectorOperations.PrintVectorToFile(contactContactivity, "Results/contactivity" + k.ToString() + ".dat");
            }
        }
    }
}
