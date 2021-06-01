using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public static class TwoBlocksInContact3D
    {
        //private const double BlockLength = 1.0;
        private const double ElementSize = 1.0;
        private static int ElementsNumber = 9;
        private const double Gap = 0.001;
        private const int nodesPerSide = 3;
        private static Dictionary<int, INode> nodes;
        public static ISolver newSolu;

        private static Dictionary<int, INode> CreateNodes()
        {

            nodes = new Dictionary<int, INode>();

            nodes[1] = new Node(0.0, 0.0, 0.0);
            nodes[2] = new Node(ElementSize, 0.0, 0.0);
            nodes[3] = new Node(2.0 * ElementSize, 0.0, 0.0);

            nodes[4] = new Node(0.0, 0.0, ElementSize);
            nodes[5] = new Node(ElementSize, 0.0, ElementSize);
            nodes[6] = new Node(2.0 * ElementSize, 0.0, ElementSize);

            nodes[7] = new Node(0.0, 0.0, 2.0 * ElementSize);
            nodes[8] = new Node(ElementSize, 0.0, 2.0 * ElementSize);
            nodes[9] = new Node(2.0 * ElementSize, 0.0, 2.0 * ElementSize);

            nodes[10] = new Node(0.0, ElementSize, 0.0);
            nodes[11] = new Node(ElementSize, ElementSize, 0.0);
            nodes[12] = new Node(2.0 * ElementSize, ElementSize, 0.0);

            nodes[13] = new Node(0.0, ElementSize, ElementSize);
            nodes[14] = new Node(ElementSize, ElementSize, ElementSize);
            nodes[15] = new Node(2.0 * ElementSize, ElementSize, ElementSize);

            nodes[16] = new Node(0.0, ElementSize, 2.0 * ElementSize);
            nodes[17] = new Node(ElementSize, ElementSize, 2.0 * ElementSize);
            nodes[18] = new Node(2.0 * ElementSize, ElementSize, 2.0 * ElementSize);

            nodes[19] = new Node(ElementSize / 2.0, ElementSize + Gap, ElementSize / 2.0);
            nodes[20] = new Node(1.5 * ElementSize, ElementSize + Gap, ElementSize / 2.0);
            nodes[21] = new Node(ElementSize / 2.0, ElementSize + Gap, 1.5 * ElementSize);
            nodes[22] = new Node(1.5 * ElementSize, ElementSize + Gap, 1.5 * ElementSize);

            nodes[23] = new Node(ElementSize / 2.0, 2.0 * ElementSize + Gap, ElementSize / 2.0);
            nodes[24] = new Node(1.5 * ElementSize, 2.0 * ElementSize + Gap, ElementSize / 2.0);
            nodes[25] = new Node(ElementSize / 2.0, 2.0 * ElementSize + Gap, 1.5 * ElementSize);
            nodes[26] = new Node(1.5 * ElementSize, 2.0 * ElementSize + Gap, 1.5 * ElementSize);

            nodes[27] = new Node(1.5 * ElementSize, 3.0 * ElementSize + Gap, ElementSize / 2.0);
            nodes[28] = new Node(ElementSize / 2.0, 3.0 * ElementSize + Gap, 1.5 * ElementSize);

            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {
            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();

            connectivity[1] = new Dictionary<int, int>() { { 1, 1 }, { 2, 2 }, { 3, 11 }, { 4, 10 }, { 5, 4 }, { 6, 5 }, { 7, 14 }, { 8, 13 } };
            connectivity[2] = new Dictionary<int, int>() { { 1, 2 }, { 2, 3 }, { 3, 12 }, { 4, 11 }, { 5, 5 }, { 6, 6 }, { 7, 15 }, { 8, 14 } };
            connectivity[3] = new Dictionary<int, int>() { { 1, 4 }, { 2, 5 }, { 3, 14 }, { 4, 13 }, { 5, 7 }, { 6, 8 }, { 7, 17 }, { 8, 16 } };
            connectivity[4] = new Dictionary<int, int>() { { 1, 5 }, { 2, 6 }, { 3, 15 }, { 4, 14 }, { 5, 8 }, { 6, 9 }, { 7, 18 }, { 8, 17 } };

            connectivity[5] = new Dictionary<int, int>() { { 1, 19 }, { 2, 20 }, { 3, 24 }, { 4, 23 }, { 5, 21 }, { 6, 22 }, { 7, 26 }, { 8, 25 } };

            connectivity[6] = new Dictionary<int, int>() { { 1, 16 }, { 2, 17 }, { 3, 14 }, { 4, 13 }, { 5, 21 } };
            connectivity[7] = new Dictionary<int, int>() { { 1, 17 }, { 2, 18 }, { 3, 15 }, { 4, 14 }, { 5, 22 } };
            connectivity[8] = new Dictionary<int, int>() { { 1, 13 }, { 2, 14 }, { 3, 11 }, { 4, 10 }, { 5, 19 } };
            connectivity[9] = new Dictionary<int, int>() { { 1, 14 }, { 2, 15 }, { 3, 12 }, { 4, 11 }, { 5, 20 } };

            connectivity[10] = new Dictionary<int, int>() { { 1, 24 }, { 2, 27 } };
            connectivity[11] = new Dictionary<int, int>() { { 1, 25 }, { 2, 28 } };

            ElementsNumber = connectivity.Count;
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= nodes.Count - 2; i++)
            {
                nodeFAT[i] = new bool[] { true, true, true, false, false, false };
            }
            nodeFAT[nodes.Count] = new bool[] { true, true, false, false, false, false };
            nodeFAT[nodes.Count - 1] = new bool[] { true, true, false, false, false, false };

            return nodeFAT;
        }

        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = 200.0 * 1e9;
            double A = 1.0;
            string type = "Hex8";
            string type2 = "ContactNtS3D";
            string type3 = "Bar2D";


            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= 5; i++)
            {
                elementProperties[i] = new ElementProperties(E, A, type);
                elementProperties[i].Density = 8000.0;
                elementProperties[i].PoissonRatio = 0.0;

                //elementProperties[i].Thickness = 0.1;
            }
            for (int i = 6; i <= 9; i++)
            {
                elementProperties[i] = new ElementProperties(E/10.0, A, type2);
            }
            for (int i = 10; i <= 11; i++)
            {
                elementProperties[i] = new ElementProperties(E/10.0, A/100.0, type3);
            }
            return elementProperties;
        }

        private static IAssembly CreateAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes();
            assembly.NodeFreedomAllocationList = CreateNodeFAT();
            assembly.ElementsConnectivity = CreateConnectivity();
            assembly.ElementsProperties = CreateElementProperties();


            assembly.BoundedDOFsVector = new int[47];
            for (int i = 1; i <= 27; i++)
            {
                assembly.BoundedDOFsVector[i - 1] = i;
            }
            int k = 1;
            for (int j = 1; j <= 8; j++)
            {
                assembly.BoundedDOFsVector[26 + k] = 54 + 3 * (j - 1) + 1;
                k += 1;
                assembly.BoundedDOFsVector[26 + k] = 54 + 3 * (j - 1) + 3;
                k += 1;
            }
            assembly.BoundedDOFsVector[43] = 79;
            assembly.BoundedDOFsVector[44] = 80;
            assembly.BoundedDOFsVector[45] = 81;
            assembly.BoundedDOFsVector[46] = 82;
            return assembly;
        }

        public static Results RunStaticExample()
        {
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();

            //ISolver newSolu = new StaticSolver();
            newSolu.LinearScheme = new LUFactorization();
            //newSolu.NonLinearScheme = new LoadControlledNewtonRaphson();
            newSolu.ActivateNonLinearSolver = true;
            newSolu.NonLinearScheme.numberOfLoadSteps = 100;

            double[] externalForces = new double[82];
            externalForces[76] = -1e8;
            //externalForces[73] = -1000.0;
            //externalForces[70] = -1000.0;
            externalForces[67] = -1e8;


            double[] reducedExternalFVector = BoundaryConditionsImposition.ReducedVector(externalForces, elementsAssembly.BoundedDOFsVector);

            newSolu.AssemblyData = elementsAssembly;
            newSolu.Solve(reducedExternalFVector);
            newSolu.PrintSolution();

            return new Results() { NonlinearSolution = new List<double[]>(), SelectedDOF = 2, SolutionType = "Nonlinear" };
        }

        public static void RunDynamicExample()
        {
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;



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