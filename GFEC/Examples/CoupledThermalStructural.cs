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
            assembly.BoundedDOFsVector = new int[] { 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165 };
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
            ShowToGUI.PlotFinalGeometry(elementsAssembly);

            double[] solVector2 = new double[280];
            List<double[]> structuralSolutions = new List<double[]>();
                int[] BoundedDOFsVector = new int[] { 1, 2, 31, 32, 61, 62, 91, 92, 121, 122, 179, 180, 209, 210, 239, 240, 269, 270, 299, 300 };
                double[] externalForces2 = new double[300];
                for (int i = 1; i <= 5; i++)
                {
                    structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson(solVector2);
                    externalForces2[135] = -10000.0 * i;
                    externalForces2[137] = -10000.0 * i;
                    externalForces2[139] = -10000.0 * i;
                    externalForces2[141] = -10000.0 * i;
                    externalForces2[143] = -10000.0 * i;
                    externalForces2[145] = -10000.0 * i;
                    externalForces2[147] = -10000.0 * i;
                    externalForces2[149] = -10000.0 * i;
                    double[] reducedExternalForces2 = BoundaryConditionsImposition.ReducedVector(externalForces2, BoundedDOFsVector);
                    structuralSolution.AssemblyData = elementsAssembly;
                    structuralSolution.Solve(reducedExternalForces2);
                    solVector2 = structuralSolution.GetSolution();
                    structuralSolutions.Add(solVector2);
                }
                Dictionary<int, double[]> intForces = structuralSolution.GetInternalForces();
                Dictionary<int, double[]> elementInternalForces = elementsAssembly.GetElementsInternalForces(structuralSolutions[0]);
                List<string> elementTypes = elementsAssembly.GetElementsType();
                double[] completeFinalSolutionVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solVector2, BoundedDOFsVector);
                Dictionary<int, INode> finalNodesList = new Dictionary<int, INode>();
                finalNodesList = Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, completeFinalSolutionVector);
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





            double[] temperatures = new double[135];
            List<double[]> thermalSolutions = new List<double[]>();
            int[] BoundedDOFsVectorForHeat = new int[] { 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165 };
            for (int i = 1; i <= 5; i++)
            {
                thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson(temperatures);
                double[] externalHeatFlux = new double[150];
                double[] reducedExternalHeatFlux = BoundaryConditionsImposition.ReducedVector(externalHeatFlux, BoundedDOFsVectorForHeat);
                thermalSolution.AssemblyData = elementsAssembly2;
                thermalSolution.Solve(reducedExternalHeatFlux);
                temperatures = thermalSolution.GetSolution();
                thermalSolutions.Add(temperatures);
            }

            List<double> X = new List<double>();
            List<double> Y = new List<double>();
            double[] Z = thermalSolutions[4];
            foreach (var node in elementsAssembly2.Nodes)
            {
                X.Add(node.Value.XCoordinate);
                Y.Add(node.Value.YCoordinate);

            }
            double[] Xvec = X.ToArray();
            double[] Yvec = Y.ToArray();

            GnuPlot.Set("terminal png size 500, 300");
            GnuPlot.Set("output 'gnuplot.png'");
            GnuPlot.Set("pm3d");
            GnuPlot.Set("dgrid3d");
            GnuPlot.Set("view map");
            GnuPlot.SPlot(Xvec, Yvec, Z);
            GnuPlot.Set("output");

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
