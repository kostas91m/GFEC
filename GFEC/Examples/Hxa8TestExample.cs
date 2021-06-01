using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public static class Hxa8TestExample
    {
        private const double BlockLength = 1.0;
        private static Dictionary<int, INode> nodes;
        public static ISolver newSolu;

        private static Dictionary<int, INode> CreateNodes()
        {
            //nodesPerSide = (int)(BlockLength / ElementSize + 1.0);
            //ElementsNumber = (int)Math.Pow(BlockLength / ElementSize, 2.0);
            nodes = new Dictionary<int, INode>();

            nodes[1] = new Node(0.0, 0.0, 0.0);
            nodes[2] = new Node(1.0, 0.0, 0.0);
            nodes[3] = new Node(1.0, 0.0, 1.0);
            nodes[4] = new Node(0.0, 0.0, 1.0);
            nodes[5] = new Node(0.0, 1.0, 0.0);
            nodes[6] = new Node(1.0, 1.0, 0.0);
            nodes[7] = new Node(1.0, 1.0, 1.0);
            nodes[8] = new Node(0.0, 1.0, 1.0);
            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {
            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();

            connectivity[1] = new Dictionary<int, int>() { { 1, 4 }, { 2, 3 }, { 3, 2 }, { 4, 1 }, { 5, 8 }, { 6, 7 }, { 7, 6 }, { 8, 5 } };
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= nodes.Count; i++)
            {
                nodeFAT[i] = new bool[] { true, true, true, false, false, false };
            }
            return nodeFAT;
        }

        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = 200.0e9;
            double A = 1.0;
            string type = "Hex8";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            elementProperties[1] = new ElementProperties(E, A, type);
            elementProperties[1].Density = 8000.0;

            return elementProperties;
        }

        private static IAssembly CreateAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes();
            assembly.ElementsConnectivity = CreateConnectivity();
            assembly.ElementsProperties = CreateElementProperties();
            assembly.NodeFreedomAllocationList = CreateNodeFAT();

            assembly.BoundedDOFsVector = new int[12];
            for (int i = 1; i <= 12; i++)
            {
                assembly.BoundedDOFsVector[i - 1] = i;
            }
            return assembly;
        }

        public static Results RunStaticExample()
        {
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();

            //ISolver newSolu = new StaticSolver();
            newSolu.LinearScheme = new CholeskyFactorization();
            //newSolu.NonLinearScheme = new LoadControlledNewtonRaphson();
            newSolu.ActivateNonLinearSolver = true;
            newSolu.NonLinearScheme.Tolerance = 1e-5;
            newSolu.NonLinearScheme.numberOfLoadSteps = 10;

            double[] externalForces = new double[24];
            externalForces[22] = -1000000000.0;


            double[] reducedExternalFVector = BoundaryConditionsImposition.ReducedVector(externalForces, elementsAssembly.BoundedDOFsVector);

            newSolu.AssemblyData = elementsAssembly;
            newSolu.Solve(reducedExternalFVector);
            double[] resultSolution = newSolu.GetSolution();
            newSolu.PrintSolution();

            return new Results() { NonlinearSolution = new List<double[]>(), SelectedDOF = 2, SolutionType = "Nonlinear" };
        }

        public static void RunDynamicExample()
        {
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = false;



            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = new double[462];
            initialValues.InitialDisplacementVector = new double[462];
            //initialValues.InitialDisplacementVector[7] = -0.02146;
            initialValues.InitialVelocityVector = new double[462];
            initialValues.InitialTime = 0.0;

            ExplicitSolver newSolver = new ExplicitSolver(1.0, 1000000);
            newSolver.Assembler = elementsAssembly;

            newSolver.InitialValues = initialValues;
            newSolver.ExternalForcesVector = new double[462];
            for (int i = 441; i <= 462; i += 2)
            {
                newSolver.ExternalForcesVector[i] = -10000.0;
            }
            newSolver.LinearSolver = new CholeskyFactorization();
            newSolver.ActivateNonLinearSolution = true;
            newSolver.SolveExplicit();
            //newSolver.PrintExplicitSolution();
        }

    }
}