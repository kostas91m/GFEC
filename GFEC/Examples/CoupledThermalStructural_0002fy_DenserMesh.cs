using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;

namespace GFEC
{
    public static class CoupledThermalStructural_0002fy_DenserMesh
    {
        private const int totalNodes = 780;
        private const int totalElements = 616;
        private const int nodesInXCoor = 78;
        private const int nodesInYCoor = 5;
        private const double scaleFactor = 1.0;
        private const double xIntervals = 0.1;
        private const double xIntervals_2 = 0.01;
        private const double yIntervals = 0.1;
        private const double offset = 0.7;
        private const double gap = 0.000000001;
        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            //Upper cantilever
            int k;
            k = 1;
            for (int j = 0; j < nodesInYCoor; j++)
            {
                for (int i = 0; i < nodesInXCoor; i++)
                {
                    if (i <= 7)
                    {
                        nodes[k] = new Node(i * xIntervals * scaleFactor, j * yIntervals * scaleFactor);
                        k += 1;
                    }
                    else
                    {
                        nodes[k] = new Node(7 * xIntervals * scaleFactor + (i - 7) * xIntervals_2 * scaleFactor, j * yIntervals * scaleFactor);
                        k += 1;
                    }

                }
            }

            //Lower cantilever
            for (int j = 0; j < nodesInYCoor; j++)
            {
                for (int i = 0; i < nodesInXCoor; i++)
                {
                    if (i <= 70)
                    {
                        nodes[k] = new Node(i * xIntervals_2 * scaleFactor + offset, j * yIntervals * scaleFactor - ((nodesInYCoor - 1) * yIntervals + gap));
                        k += 1;
                    }
                    else
                    {
                        nodes[k] = new Node((i - 70) * xIntervals * scaleFactor + offset + 70 * xIntervals_2, j * yIntervals * scaleFactor - ((nodesInYCoor - 1) * yIntervals + gap));
                        k += 1;
                    }

                }
            }
            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int k = 1;
            for (int j = 0; j <= nodesInYCoor - 2; j++)
            {
                for (int i = 1; i <= nodesInXCoor - 1; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + j * nodesInXCoor }, { 2, i + 1 + j * nodesInXCoor }, { 3, i + 1 + nodesInXCoor + j * nodesInXCoor }, { 4, i + nodesInXCoor + j * nodesInXCoor } };
                    k += 1;
                }

            }

            for (int j = 0 + nodesInYCoor; j <= nodesInYCoor - 2 + nodesInYCoor; j++)
            {
                for (int i = 1; i <= nodesInXCoor - 1; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + j * nodesInXCoor }, { 2, i + 1 + j * nodesInXCoor }, { 3, i + 1 + nodesInXCoor + j * nodesInXCoor }, { 4, i + nodesInXCoor + j * nodesInXCoor } };
                    k += 1;
                }

            }

            //Contact elements
            //int xContactNodesNumber = (int)(nodesInXCoor - offset / 0.1);
            //List<int> lowerContactNodes = new List<int>();
            //List<int> upperContactNodes = new List<int>();

            //for (int i = nodesInXCoor - xContactNodesNumber + 1; i <= nodesInXCoor; i++)
            //{

            //}
            int n = 8;
            int m = 703;
            for (int i = 617; i <= 687; i++)
            {
                connectivity[i] = new Dictionary<int, int>() { { 1, m }, { 2, n } };
                m += 1;
                n += 1;
            }
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= totalNodes; i++)
            {
                nodeFAT[i] = new bool[] { true, true, false, false, false, false };
            }
            return nodeFAT;
        }

        private static Dictionary<int, bool[]> CreateThermalNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= totalNodes; i++)
            {
                nodeFAT[i] = new bool[] { true, false, false, false, false, false };
            }
            return nodeFAT;
        }

        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = 200.0e9;
            double A = 0.01;
            string type = "Quad4";
            string type2 = "ContactNtN2D";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= totalElements; i++)
            {
                elementProperties[i] = new ElementProperties(E, A, type);
            }

            for (int i = 1; i <= totalElements; i++)
            {
                elementProperties[i].Density = 7850.0;
                elementProperties[i].Thickness = 0.1;
            }

            for (int i = 617; i <= 687; i++)
            {
                elementProperties[i] = new ElementProperties(E / 1000.0, A, type2);
                elementProperties[i].Density = 7850.0;
                elementProperties[i].Thickness = 0.1;
            }
            return elementProperties;
        }

        private static Dictionary<int, IElementProperties> CreateThermalElementProperties()
        {
            double thermalCond = 15.5;
            double A = 0.001;
            double b = 0.1;
            double a1 = 0.1;
            double a2 = 0.01;
            string type = "Quad4Th2";
            string type2 = "ContactNtN2DTh";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= totalElements; i++)
            {
                if (i <= 7 | (i > 77 && i <= 84) | (i > 154 && i <= 161) | (i > 231 && i <= 238) | (i > 308 && i <= 315) | (i > 385 && i <= 392) | (i > 462 && i <= 469) | (i > 539 && i <= 546))
                {
                    elementProperties[i] = new ElementProperties();
                    elementProperties[i].ElementType = type;
                    elementProperties[i].ThermalConductivity = thermalCond;
                    elementProperties[i].lx = a1;
                    elementProperties[i].ly = b;
                }
                else
                {
                    elementProperties[i] = new ElementProperties();
                    elementProperties[i].ElementType = type;
                    elementProperties[i].ThermalConductivity = thermalCond;
                    elementProperties[i].lx = a2;
                    elementProperties[i].ly = b;
                }

            }
            for (int i = 617; i <= 687; i++)
            {
                elementProperties[i] = new ElementProperties();
                elementProperties[i].ElementType = type2;
                elementProperties[i].SectionArea = A;
                elementProperties[i].ContactForceValue = 0.0;
            }
            elementProperties[617].SectionArea = 0.011;
            elementProperties[687].SectionArea = 0.0005;
            return elementProperties;
        }

        private static IAssembly CreateAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes();
            assembly.ElementsConnectivity = CreateConnectivity();
            assembly.ElementsProperties = CreateElementProperties();
            assembly.NodeFreedomAllocationList = CreateNodeFAT();
            assembly.BoundedDOFsVector = new int[] { 1, 2, 157, 158, 313, 314, 469, 470, 625, 626, 782, 935, 936, 1091, 1092, 1247, 1248, 1403, 1404, 1559, 1560 };
            return assembly;
        }

        private static IAssembly CreateThermalAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes();
            assembly.ElementsConnectivity = CreateConnectivity();
            assembly.ElementsProperties = CreateThermalElementProperties();
            assembly.NodeFreedomAllocationList = CreateThermalNodeFAT();
            assembly.BoundedDOFsVector = new int[] { 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468 };
            return assembly;
        }

        public static Results RunStaticExample()
        {
            #region Structural
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();

            //Gnuplot graphs
            ShowToGUI.PlotInitialGeometry(elementsAssembly);


            ISolver structuralSolution = new StaticSolver();
            structuralSolution.LinearScheme = new LUFactorization();
            structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 10;
            //int[] BoundedDOFsVector2 = new int[] { 1, 2, 31, 32, 61, 62, 91, 92, 121, 122, 179, 180, 209, 210, 239, 240, 269, 270, 299, 300 };
            double[] externalForces3 = new double[1560];
            externalForces3[627] = -2500.0;
            externalForces3[629] = -5000.0;
            externalForces3[631] = -5000.0;
            externalForces3[633] = -5000.0;
            externalForces3[635] = -5000.0;
            externalForces3[637] = -5000.0;
            externalForces3[639] = -2750.0;
            externalForces3[779] = -250.0;
            for (int i = 641; i <= 777; i += 2)
            {
                externalForces3[i] = -500.0;

            }
            double[] reducedExternalForces3 = BoundaryConditionsImposition.ReducedVector(externalForces3, elementsAssembly.BoundedDOFsVector);
            structuralSolution.AssemblyData = elementsAssembly;
            structuralSolution.Solve(reducedExternalForces3);
            double[] solvector3 = structuralSolution.GetSolution();
            elementsAssembly.UpdateDisplacements(solvector3);
            ShowToGUI.PlotFinalGeometry(elementsAssembly);
            double[] fullSolVector3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector3, elementsAssembly.BoundedDOFsVector);
            Dictionary<int, INode> finalNodes = Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullSolVector3);
            double[] xFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item1;
            double[] yFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item2;
            Dictionary<int, double[]> allStepsSolutions = structuralSolution.GetAllStepsSolutions();

            Dictionary<int, Dictionary<int, double[]>> allStepsContactForces = new Dictionary<int, Dictionary<int, double[]>>();
            Dictionary<int, double[]> elementsInternalContactForcesVector;
            for (int i = 1; i <= allStepsSolutions.Count; i++)
            {
                elementsInternalContactForcesVector = new Dictionary<int, double[]>();
                elementsAssembly.UpdateDisplacements(allStepsSolutions[i]);
                for (int j = 617; j <= 687; j++)
                {
                    elementsInternalContactForcesVector[j] = elementsAssembly.ElementsAssembly[j].CreateInternalGlobalForcesVector();
                }
                allStepsContactForces[i] = elementsInternalContactForcesVector;
            }



            //    double[] solVector2 = new double[280];
            List<double[]> structuralSolutions = new List<double[]>();
            //    int[] BoundedDOFsVector = new int[] { 1, 2, 31, 32, 61, 62, 91, 92, 121, 122, 179, 180, 209, 210, 239, 240, 269, 270, 299, 300 };
            //    double[] externalForces2 = new double[300];
            //    for (int i = 1; i <= 5; i++)
            //    {
            //        structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson(solVector2);
            //        externalForces2[135] = -10000.0 * i;
            //        externalForces2[137] = -10000.0 * i;
            //        externalForces2[139] = -10000.0 * i;
            //        externalForces2[141] = -10000.0 * i;
            //        externalForces2[143] = -10000.0 * i;
            //        externalForces2[145] = -10000.0 * i;
            //        externalForces2[147] = -10000.0 * i;
            //        externalForces2[149] = -10000.0 * i;
            //        double[] reducedExternalForces2 = BoundaryConditionsImposition.ReducedVector(externalForces2, BoundedDOFsVector);
            //        structuralSolution.AssemblyData = elementsAssembly;
            //        structuralSolution.Solve(reducedExternalForces2);
            //        solVector2 = structuralSolution.GetSolution();
            //        structuralSolutions.Add(solVector2);
            //    }
            //    Dictionary<int, double[]> intForces = structuralSolution.GetInternalForces();
            //    Dictionary<int, double[]> elementInternalForces = elementsAssembly.GetElementsInternalForces(structuralSolutions[0]);
            //    List<string> elementTypes = elementsAssembly.GetElementsType();
            //    double[] completeFinalSolutionVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solVector2, BoundedDOFsVector);
            //    Dictionary<int, INode> finalNodesList = new Dictionary<int, INode>();
            //    finalNodesList = Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, completeFinalSolutionVector);
            #endregion


            #region Thermal
            List<double[]> thermalSolutions = new List<double[]>();
            for (int k = 1; k <= allStepsSolutions.Count; k++)
            {
                IAssembly elementsAssembly2 = CreateThermalAssembly();

                for (int j = 617; j < 687; j++)
                {
                    double[] contactForce = allStepsContactForces[k][j];
                    elementsAssembly2.ElementsProperties[j].ContactForceValue = VectorOperations.VectorNorm2(new double[] { contactForce[2], contactForce[3] });
                }

                elementsAssembly2.CreateElementsAssembly();
                elementsAssembly2.ActivateBoundaryConditions = true;
                double[,] globalStiffnessMatrix2 = elementsAssembly2.CreateTotalStiffnessMatrix();

                ISolver thermalSolution = new StaticSolver();
                thermalSolution.LinearScheme = new CholeskyFactorization();
                thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                thermalSolution.ActivateNonLinearSolver = true;
                thermalSolution.NonLinearScheme.numberOfLoadSteps = 10;

                thermalSolution.AssemblyData = elementsAssembly2;
                double[] externalHeatFlux = new double[780];
                for (int i = 312; i <= 389; i++)
                {
                    if (i > 312 && i <= 318)
                    {
                        externalHeatFlux[i] = 100.0;
                    }
                    else
                    {
                        externalHeatFlux[i] = 10.0;
                    }
                    externalHeatFlux[312] = 50.0;
                    externalHeatFlux[319] = 55.0;
                    externalHeatFlux[389] = 5.0;

                }
                double[] reducedExternalHeatFlux = BoundaryConditionsImposition.ReducedVector(externalHeatFlux, thermalSolution.AssemblyData.BoundedDOFsVector);
                thermalSolution.Solve(reducedExternalHeatFlux);
                double[] tempSol = thermalSolution.GetSolution();
                thermalSolutions.Add(tempSol);
            }

            int[] thermalBoundCond = new int[] { 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468 };
            double[] fullStructuralSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[2], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[4], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[6], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[8], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[10], elementsAssembly.BoundedDOFsVector);
            double[] fullThermalSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[1], thermalBoundCond);
            double[] fullThermalSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[3], thermalBoundCond);
            double[] fullThermalSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[5], thermalBoundCond);
            double[] fullThermalSol4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[7], thermalBoundCond);
            double[] fullThermalSol5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[9], thermalBoundCond);
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol1), fullThermalSol1, @"C:\Users\Public\Documents\Results1_0.002fy_DM.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol2), fullThermalSol2, @"C:\Users\Public\Documents\Results2_0.002fy_DM.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol3), fullThermalSol3, @"C:\Users\Public\Documents\Results3_0.002fy_DM.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol4), fullThermalSol4, @"C:\Users\Public\Documents\Results4_0.002fy_DM.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol5), fullThermalSol5, @"C:\Users\Public\Documents\Results5_0.002fy_DM.dat");



            //double[] temperatures = new double[135];

            //int[] BoundedDOFsVectorForHeat = new int[] { 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165 };
            //for (int i = 1; i <= 5; i++)
            //{
            //    thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson(temperatures);
            //    double[] externalHeatFlux = new double[150];
            //    double[] reducedExternalHeatFlux = BoundaryConditionsImposition.ReducedVector(externalHeatFlux, BoundedDOFsVectorForHeat);
            //    thermalSolution.AssemblyData = elementsAssembly2;
            //    thermalSolution.Solve(reducedExternalHeatFlux);
            //    temperatures = thermalSolution.GetSolution();
            //    thermalSolutions.Add(temperatures);
            //}
            //List<double> X = new List<double>();
            //List<double> Y = new List<double>();
            //double[] fullTempSol = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(tempSol, thermalSolution.AssemblyData.BoundedDOFsVector);
            //double[] Z = fullTempSol;//thermalSolutions[4];
            //foreach (var node in elementsAssembly2.Nodes)
            //{
            //    X.Add(node.Value.XCoordinate);
            //    Y.Add(node.Value.YCoordinate);

            //}
            //double[] Xvec = X.ToArray();
            //double[] Yvec = Y.ToArray();


            //double[] Xvec1 = new double[75];
            //double[] Yvec1 = new double[75];
            //double[] Zvec1 = new double[75];
            //double[] Xvec2 = new double[75];
            //double[] Yvec2 = new double[75];
            //double[] Zvec2 = new double[75];
            //Array.Copy(Xvec, 75, Xvec2, 0, 75);
            //Array.Copy(Yvec, 75, Yvec2, 0, 75);
            //Array.Copy(Z, 75, Zvec2, 0, 75);
            //Array.Copy(Xvec, 0, Xvec1, 0, 75);
            //Array.Copy(Yvec, 0, Yvec1, 0, 75);
            //Array.Copy(Z, 0, Zvec1, 0, 75);
            //GnuPlot.Set("terminal png size 500, 300");
            //GnuPlot.Set("output 'gnuplot.png'");

            //GnuPlot.HoldOn();
            //GnuPlot.Set("cbrange[0:20.0]");
            //GnuPlot.Set("palette defined(0 \"blue\", 1 \"red\")");
            //GnuPlot.Set("pm3d");
            //GnuPlot.Set("dgrid3d");
            //GnuPlot.Set("view map");
            //GnuPlot.SPlot(Xvec1, Yvec1, Zvec1);
            //GnuPlot.SPlot(Xvec2, Yvec2, Zvec2);
            //GnuPlot.Set("output");
            //List<HeatMapData> plots = new List<HeatMapData>();
            //plots.Add(new HeatMapData() { Xcoordinates = Xvec1, Ycoordinates = Yvec1, Temperatures = Zvec1 });
            //plots.Add(new HeatMapData() { Xcoordinates = Xvec2, Ycoordinates = Yvec2, Temperatures = Zvec2 });
            ////ShowToGUI.PlotHeatMap(plots);

            //double[] Xvec1Final = new double[75];
            //double[] Yvec1Final = new double[75];
            //double[] Xvec2Final = new double[75];
            //double[] Yvec2Final = new double[75];

            //Array.Copy(xFinalNodalCoor, 0, Xvec1Final, 0, 75);
            //Array.Copy(yFinalNodalCoor, 0, Yvec1Final, 0, 75);
            //Array.Copy(xFinalNodalCoor, 75, Xvec2Final, 0, 75);
            //Array.Copy(yFinalNodalCoor, 75, Yvec2Final, 0, 75);

            //List<HeatMapData> plots2 = new List<HeatMapData>();
            //plots2.Add(new HeatMapData() { Xcoordinates = Xvec1Final, Ycoordinates = Yvec1Final, Temperatures = Zvec1 });
            //plots2.Add(new HeatMapData() { Xcoordinates = Xvec2Final, Ycoordinates = Yvec2Final, Temperatures = Zvec2 });
            //GnuPlot.HoldOn();
            //GnuPlot.Set("pm3d");
            //GnuPlot.Set("dgrid3d");
            //GnuPlot.Set("view map");
            //GnuPlot.SPlot(new double[] { -1.0, 2.0, 1.0, -1.0}, new double[] { 1.0, 2.0, -1.0, 1.0 }, new double[] { 2, 1, 3, 2 });
            ////GnuPlot.SPlot(new double[] { -1.0, 1.0, 3.0 }, new double[] { 2.0, 2.0, -1.0 }, new double[] { 5, 4, 9 });
            ////GnuPlot.Plot(Xvec2Final, Yvec2Final);
            //ShowToGUI.PlotHeatMap(plots2);

            //ExportToFile.ExportGeometryDataWithTemperatures(finalNodes, fullTempSol);


            GnuPlot.Close();

            while (true)
            {
                if (File.Exists(AppContext.BaseDirectory + "gnuplot.png") && new FileInfo(AppContext.BaseDirectory + "gnuplot.png").Length > 0)
                {
                    break;

                }
                Thread.Sleep(100);
            }
            GnuPlot.KillProcess();
            #endregion
            return new Results() { NonlinearSolution = structuralSolutions, SelectedDOF = 2, SolutionType = "Nonlinear" };
        }

        public static void RunDynamicExample()
        {
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;

            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = new double[6];
            initialValues.InitialDisplacementVector = new double[6];
            //initialValues.InitialDisplacementVector[7] = -0.02146;
            initialValues.InitialVelocityVector = new double[6];
            initialValues.InitialTime = 0.0;

            ExplicitSolver newSolver = new ExplicitSolver(1.0, 10000);
            newSolver.Assembler = elementsAssembly;

            newSolver.InitialValues = initialValues;
            newSolver.ExternalForcesVector = new double[] { 0, 0, 0, 0, -50000, -50000 };
            newSolver.LinearSolver = new CholeskyFactorization();
            newSolver.ActivateNonLinearSolution = true;
            newSolver.SolveNewmark();
            newSolver.PrintExplicitSolution();//
        }

    }
}
