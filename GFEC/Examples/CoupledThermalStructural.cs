using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public static class CoupledThermalStructural
    {
        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            double[] X = new double[100];
            double[] Y = new double[100];
            double[] Z = new double[100];
            double scaleFactor = 1.0;
            double xIntervals = 0.1;
            double yIntervals = 0.1;
            int k, m;
            m = 1;
            k = 1;
            for (int j = 0; j < 11; j++)
            {
                //Y[m - 1] = j * yIntervals * scaleFactor;
                for (int i = 0; i < 11; i++)
                {
                    nodes[k] = new Node(i * xIntervals * scaleFactor, j * yIntervals * scaleFactor);
                    //X[k - 1] = i * xIntervals * scaleFactor;
                    k += 1;
                }
                m += 1;
            }
            

            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {
            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int k = 1;
            for (int j = 0; j <= 9; j++)
            {
                for (int i = 1; i <= 10; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i+j*11 }, { 2, i + 1 +j*11}, { 3, i + 1 + 11+j*11 }, { 4, i + 11+j*11 } };
                    k += 1;
                }
                
            }

            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= 121; i++)
            {
                nodeFAT[i] = new bool[] { true, true, false, false, false, false };
            }
            return nodeFAT;
        }

        private static Dictionary<int, bool[]> CreateThermalNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= 121; i++)
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
            for (int i = 1; i <= 100; i++)
            {
                elementProperties[i] = new ElementProperties(E, A, type);
            }
            //elementProperties[1] = new ElementProperties(E, A, type);
            //elementProperties[2] = new ElementProperties(E, A, type);
            //elementProperties[3] = new ElementProperties(E / 1000.0, A, type2);
            //elementProperties[4] = new ElementProperties(E / 1000.0, A, type2);
            for (int i = 1; i <= 2; i++)
            {
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
            for (int i = 1; i <= 100; i++)
            {
                elementProperties[i] = new ElementProperties();
            }
            
            elementProperties[1].ElementType = type;
            elementProperties[2].ElementType = type;
            elementProperties[3].ElementType = type2;
            elementProperties[4].ElementType = type2;
            elementProperties[1].ThermalConductivity = thermalCond;
            elementProperties[2].ThermalConductivity = thermalCond;
            elementProperties[3].SectionArea = A;
            elementProperties[4].SectionArea = A;

            return elementProperties;
        }

        private static IAssembly CreateAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes();
            assembly.ElementsConnectivity = CreateConnectivity();
            assembly.ElementsProperties = CreateElementProperties();
            assembly.NodeFreedomAllocationList = CreateNodeFAT();
            assembly.BoundedDOFsVector = new int[] { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15 };
            return assembly;
        }

        private static IAssembly CreateThermalAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes();
            assembly.ElementsConnectivity = CreateConnectivity();
            assembly.ElementsProperties = CreateThermalElementProperties();
            assembly.NodeFreedomAllocationList = CreateThermalNodeFAT();
            assembly.BoundedDOFsVector = new int[] { 1, 2 };
            return assembly;
        }

        public static Results RunStaticExample()
        {
            ///Structural settings
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();

            //Gnuplot graphs
            ShowToGUI.PlotWithGnuPlot(elementsAssembly);


            ISolver newSolu = new StaticSolver();
            newSolu.LinearScheme = new PCGSolver();
            newSolu.NonLinearScheme = new LoadControlledNewtonRaphson();
            newSolu.ActivateNonLinearSolver = true;
            newSolu.NonLinearScheme.numberOfLoadSteps = 1;

            //Thermal Settings
            IAssembly elementsAssembly2 = CreateThermalAssembly();
            elementsAssembly2.CreateElementsAssembly();
            elementsAssembly2.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix2 = elementsAssembly2.CreateTotalStiffnessMatrix();

            ISolver thermalSolution = new StaticSolver();
            thermalSolution.LinearScheme = new CholeskyFactorization();
            thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
            thermalSolution.ActivateNonLinearSolver = true;
            thermalSolution.NonLinearScheme.numberOfLoadSteps = 10;

            //double[] externalFlux = new double[] { 0, 0, 0, 0, 250.0, 250.0 };
            //newSolu2.AssemblyData = elementsAssembly;
            //newSolu2.Solve(externalForces);
            //newSolu2.PrintSolution();
            //double[] tempSolution = newSolu.GetSolution();

            double[] solVector2 = new double[6];
            List<double[]> structuralSolutions = new List<double[]>();
            for (int i = 1; i <= 5; i++)
            {
                newSolu.NonLinearScheme = new LoadControlledNewtonRaphson(solVector2);
                double[] externalForces2 = new double[] { 0, 0, 0, 0, -10000.0 * i, -10000.0 * i };
                newSolu.AssemblyData = elementsAssembly;
                newSolu.Solve(externalForces2);
                solVector2 = newSolu.GetSolution();
                structuralSolutions.Add(solVector2);
            }
            Dictionary<int, double[]> intForces = newSolu.GetInternalForces();
            Dictionary<int, double[]> elementInternalForces = elementsAssembly.GetElementsInternalForces(structuralSolutions[0]);
            List<string> elementTypes = elementsAssembly.GetElementsType();

            double[] temperatures = new double[6];
            List<double[]> thermalSolutions = new List<double[]>();
            for (int i = 1; i <= 5; i++)
            {
                thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson(temperatures);
                double[] externalHeatFlux = new double[] { 0, 0, 0, 0, 250.0, 250.0 };
                thermalSolution.AssemblyData = elementsAssembly2;
                thermalSolution.Solve(externalHeatFlux);
                temperatures = thermalSolution.GetSolution();
                thermalSolutions.Add(temperatures);
            }


            //double[] completeFinalSolutionVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solVector2, new int[] { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15 });            Dictionary<int, INode> finalNodesList = new Dictionary<int, INode>();
            //finalNodesList = Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, completeFinalSolutionVector);

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
