using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;

namespace GFEC
{
    public static class CoupledThermalStructural
    {
        private const int totalNodes = 150;
        private const int totalElements = 112;
        private const int nodesInXCoor = 15;
        private const int nodesInYCoor = 5;
        private const double scaleFactor = 1.0;
        private const double xIntervals = 0.1;
        private const double yIntervals = 0.1;
        private const double offset = 0.7;
        private const double gap = 0.01;

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
                    nodes[k] = new Node(i * xIntervals * scaleFactor, j * yIntervals * scaleFactor);
                    k += 1;
                }
            }

            //Lower cantilever
            for (int j = 0; j < nodesInYCoor; j++)
            {
                for (int i = 0; i < nodesInXCoor; i++)
                {
                    nodes[k] = new Node(i * xIntervals * scaleFactor + offset, j * yIntervals * scaleFactor - ((nodesInYCoor - 1) * yIntervals + gap));
                    k += 1;
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
            connectivity[113] = new Dictionary<int, int>() { { 1, 136 }, { 2, 8 } };
            connectivity[114] = new Dictionary<int, int>() { { 1, 137 }, { 2, 9 } };
            connectivity[115] = new Dictionary<int, int>() { { 1, 138 }, { 2, 10 } };
            connectivity[116] = new Dictionary<int, int>() { { 1, 139 }, { 2, 11 } };
            connectivity[117] = new Dictionary<int, int>() { { 1, 140 }, { 2, 12 } };
            connectivity[118] = new Dictionary<int, int>() { { 1, 141 }, { 2, 13 } };
            connectivity[119] = new Dictionary<int, int>() { { 1, 142 }, { 2, 14 } };
            connectivity[120] = new Dictionary<int, int>() { { 1, 143 }, { 2, 15 } };
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
                elementProperties[i].Density = 8000.0;
                elementProperties[i].Thickness = 0.1;
            }

            for (int i = 113; i <= 120; i++)
            {
                elementProperties[i] = new ElementProperties(E, A, type2);
                elementProperties[i].Density = 8000.0;
                elementProperties[i].Thickness = 0.1;
            }
            return elementProperties;
        }

        private static Dictionary<int, IElementProperties> CreateThermalElementProperties()
        {
            double thermalCond = 60.5;
            double A = 0.5;
            string type = "Quad4Th";
            string type2 = "ContactNtN2DTh";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= totalElements; i++)
            {
                elementProperties[i] = new ElementProperties();
                elementProperties[i].ElementType = type;
                elementProperties[i].ThermalConductivity = thermalCond;
            }
            for (int i = 113; i <= 120; i++)
            {
                elementProperties[i] = new ElementProperties();
                elementProperties[i].ElementType = type2;
                elementProperties[i].ThermalConductivity = thermalCond;
                elementProperties[i].SectionArea = A;
            }
            return elementProperties;
        }

        private static IAssembly CreateAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes();
            assembly.ElementsConnectivity = CreateConnectivity();
            assembly.ElementsProperties = CreateElementProperties();
            assembly.NodeFreedomAllocationList = CreateNodeFAT();
            assembly.BoundedDOFsVector = new int[] { 1, 2, 31, 32, 61, 62, 91, 92, 121, 122, 179, 180, 209, 210, 239, 240, 269, 270, 299, 300 };
            return assembly;
        }

        private static IAssembly CreateThermalAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes();
            assembly.ElementsConnectivity = CreateConnectivity();
            assembly.ElementsProperties = CreateThermalElementProperties();
            assembly.NodeFreedomAllocationList = CreateThermalNodeFAT();
            assembly.BoundedDOFsVector = new int[] { 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90 };
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
            //ShowToGUI.PlotInitialGeometry(elementsAssembly);


            ISolver structuralSolution = new StaticSolver();
            structuralSolution.LinearScheme = new LUFactorization();
            structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 10;
            int[] BoundedDOFsVector2 = new int[] { 1, 2, 31, 32, 61, 62, 91, 92, 121, 122, 179, 180, 209, 210, 239, 240, 269, 270, 299, 300 };
            double[] externalForces3 = new double[300];
            externalForces3[135] = -100000000.0;
            externalForces3[137] = -100000000.0;
            externalForces3[139] = -100000000.0;
            externalForces3[141] = -100000000.0;
            externalForces3[143] = -100000000.0;
            externalForces3[145] = -100000000.0;
            externalForces3[147] = -100000000.0;
            externalForces3[149] = -100000000.0;
            double[] reducedExternalForces3 = BoundaryConditionsImposition.ReducedVector(externalForces3, BoundedDOFsVector2);
            structuralSolution.AssemblyData = elementsAssembly;
            structuralSolution.Solve(reducedExternalForces3);
            double[] solvector3 = structuralSolution.GetSolution();
            elementsAssembly.UpdateDisplacements(solvector3);
            //ShowToGUI.PlotFinalGeometry(elementsAssembly);
            double[] fullSolVector3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector3, elementsAssembly.BoundedDOFsVector);
            Dictionary<int, INode> finalNodes = Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullSolVector3);
            double[] xFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item1;
            double[] yFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item2;

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
            IAssembly elementsAssembly2 = CreateThermalAssembly();
            elementsAssembly2.CreateElementsAssembly();
            elementsAssembly2.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix2 = elementsAssembly2.CreateTotalStiffnessMatrix();

            ISolver thermalSolution = new StaticSolver();
            thermalSolution.LinearScheme = new CholeskyFactorization();
            thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
            thermalSolution.ActivateNonLinearSolver = true;
            thermalSolution.NonLinearScheme.numberOfLoadSteps = 10;

            thermalSolution.AssemblyData = elementsAssembly2;
            double[] externalHeatFlux = new double[150];
            for (int i = 61; i <= 75; i++)
            {
                externalHeatFlux[61] = 250.0;
            }
            double[] reducedExternalHeatFlux = BoundaryConditionsImposition.ReducedVector(externalHeatFlux, thermalSolution.AssemblyData.BoundedDOFsVector);
            thermalSolution.Solve(reducedExternalHeatFlux);
            double[] tempSol = thermalSolution.GetSolution();



            //double[] temperatures = new double[135];
            List<double[]> thermalSolutions = new List<double[]>();
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

            List<double> X = new List<double>();
            List<double> Y = new List<double>();
            double[] fullTempSol = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(tempSol, thermalSolution.AssemblyData.BoundedDOFsVector);
            double[] Z = fullTempSol;//thermalSolutions[4];
            foreach (var node in elementsAssembly2.Nodes)
            {
                X.Add(node.Value.XCoordinate);
                Y.Add(node.Value.YCoordinate);

            }
            double[] Xvec = X.ToArray();
            double[] Yvec = Y.ToArray();


            double[] Xvec1 = new double[75];
            double[] Yvec1 = new double[75];
            double[] Zvec1 = new double[75];
            double[] Xvec2 = new double[75];
            double[] Yvec2 = new double[75];
            double[] Zvec2 = new double[75];
            Array.Copy(Xvec, 75, Xvec2, 0, 75);
            Array.Copy(Yvec, 75, Yvec2, 0, 75);
            Array.Copy(Z, 75, Zvec2, 0, 75);
            Array.Copy(Xvec, 0, Xvec1, 0, 75);
            Array.Copy(Yvec, 0, Yvec1, 0, 75);
            Array.Copy(Z, 0, Zvec1, 0, 75);
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
            List<HeatMapData> plots = new List<HeatMapData>();
            plots.Add(new HeatMapData() { Xcoordinates = Xvec1, Ycoordinates = Yvec1, Temperatures = Zvec1 });
            plots.Add(new HeatMapData() { Xcoordinates = Xvec2, Ycoordinates = Yvec2, Temperatures = Zvec2 });
            //ShowToGUI.PlotHeatMap(plots);

            double[] Xvec1Final = new double[75];
            double[] Yvec1Final = new double[75];
            double[] Xvec2Final = new double[75];
            double[] Yvec2Final = new double[75];
            
            Array.Copy(xFinalNodalCoor, 0, Xvec1Final, 0, 75);
            Array.Copy(yFinalNodalCoor, 0, Yvec1Final, 0, 75);
            Array.Copy(xFinalNodalCoor, 75, Xvec2Final, 0, 75);
            Array.Copy(yFinalNodalCoor, 75, Yvec2Final, 0, 75);

            List<HeatMapData> plots2 = new List<HeatMapData>();
            plots2.Add(new HeatMapData() { Xcoordinates = Xvec1Final, Ycoordinates = Yvec1Final, Temperatures = Zvec1 });
            plots2.Add(new HeatMapData() { Xcoordinates = Xvec2Final, Ycoordinates = Yvec2Final, Temperatures = Zvec2 });
            GnuPlot.HoldOn();
            GnuPlot.Set("pm3d");
            GnuPlot.Set("dgrid3d");
            GnuPlot.Set("view map");
            GnuPlot.SPlot(new double[] { -1.0, 2.0, 1.0, -1.0}, new double[] { 1.0, 2.0, -1.0, 1.0 }, new double[] { 2, 1, 3, 2 });
            //GnuPlot.SPlot(new double[] { -1.0, 1.0, 3.0 }, new double[] { 2.0, 2.0, -1.0 }, new double[] { 5, 4, 9 });
            //GnuPlot.Plot(Xvec2Final, Yvec2Final);
            ShowToGUI.PlotHeatMap(plots2);

            double kati = 1;
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
