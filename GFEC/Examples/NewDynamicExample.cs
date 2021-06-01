using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class NewDynamicExample
    {
        private const double density = 1500;
        private const double density2 = 800;
        private const double YoungMod1 = 30 * 1e9;
        private const double YoungMod2 = 200 * 1e9;
        private const double area = 0.1;
        private const double thickness = 0.1;
        private const double xIntervals = 0.5;
        private const double yIntervals = 0.5;
        private const double barLenght = 0.5;
        private const double externalForce = -100;
        static int[] structuralBoundaryConditions;
        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;
        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            nodes[1] = new Node(0, 0);
            nodes[2] = new Node(xIntervals, 0);
            nodes[3] = new Node(2 * xIntervals, 0);
            nodes[4] = new Node(0, yIntervals);
            nodes[5] = new Node(xIntervals, yIntervals);
            nodes[6] = new Node(2 * xIntervals, yIntervals);
            nodes[7] = new Node(0, 2 * yIntervals);
            nodes[8] = new Node(xIntervals, 2 * yIntervals);
            nodes[9] = new Node(2 * xIntervals, 2 * yIntervals);
            nodes[10] = new Node(0, 2 * yIntervals + barLenght);
            nodes[11] = new Node(xIntervals, 2 * yIntervals + barLenght);
            nodes[12] = new Node(2 * xIntervals, 2 * yIntervals + barLenght);
            return nodes;
        }
        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= 12; i++)
            {
                nodeFAT[i] = new bool[] { true, true, false, false, false, false };
            }
            return nodeFAT;
        }
        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            for (int i = 1; i <= 7; i += 3)
            {
                boundedDofs.Add(2 * i - 1); //left side support
            }
            for (int i = 3; i <= 9; i += 3)
            {
                boundedDofs.Add(2 * i - 1); //right side support
            }
            for (int i = 10; i <= 12; i++)
            {
                boundedDofs.Add(2 * i - 1); //bar elements supports
                boundedDofs.Add(2 * i); //bar elements supports
            }
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }
        private static double[] CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            for (int i = 1; i <= 3; i++)
            {
                loadedStructuralDOFs.Add(i * 2);
            }
            externalForcesStructuralVector = new double[12 * 2];
            return externalForcesStructuralVector;
        }
        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            connectivity[1] = new Dictionary<int, int>() { { 1, 1 }, { 2, 2 }, { 3, 5 }, { 4, 4} };
            connectivity[2] = new Dictionary<int, int>() { { 1, 2 }, { 2, 3 }, { 3, 6 }, { 4, 5 } };
            connectivity[3] = new Dictionary<int, int>() { { 1, 4 }, { 2, 5 }, { 3, 8 }, { 4, 7 } };
            connectivity[4] = new Dictionary<int, int>() { { 1, 5 }, { 2, 6 }, { 3, 9 }, { 4, 8 } };
            connectivity[5] = new Dictionary<int, int>() { { 1, 7 }, { 2, 10 }};
            connectivity[6] = new Dictionary<int, int>() { { 1, 8 }, { 2, 11 } };
            connectivity[7] = new Dictionary<int, int>() { { 1, 9 }, { 2, 12 } };

            return connectivity;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E1 = YoungMod1;
            double E2 = YoungMod2;
            double A = area;
            string type = "Quad4";
            string type2 = "Bar2D";
            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= 4; i++)
            {
                elementProperties[i] = new ElementProperties(E1, A, type);
            }
            for (int i = 1; i <= 4; i++)
            {
                elementProperties[i].Density = density;
                elementProperties[i].Thickness = thickness;
            }
            for (int i = 5; i <= 7; i++)
            {
                elementProperties[i] = new ElementProperties(E2, A, type2);
            }
            for (int i = 5; i <= 7; i++)
            {
                elementProperties[i].Density = density2;
                elementProperties[i].Thickness = thickness;
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
            CreateStructuralBoundaryConditions();
            CreateStructuralLoadVector();
            assembly.BoundedDOFsVector = structuralBoundaryConditions;
            return assembly;
        }
        public static Results RunExample()
        {
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;
            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = new double[] { 10, 0.0, 10, 10, 10, 0.0, 10, 10, 10, 0.0, 10, 10 };
            initialValues.InitialDisplacementVector = new double[] { -0.2, 0.0, -0.2, -0.2, -0.2, 0.0, -0.2, -0.2, -0.2, 0.0, -0.2, -0.2 };
            initialValues.InitialVelocityVector = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            initialValues.InitialTime = 0.0;

            ExplicitSolver newSolver = new ExplicitSolver(10, 100);
            newSolver.Assembler = elementsAssembly;
            double[] externalForces = externalForcesStructuralVector;
            foreach (var dof in loadedStructuralDOFs)
            {
                externalForces[dof - 1] = externalForce;

            }
            newSolver.InitialValues = initialValues;
            newSolver.ExternalForcesVector = BoundaryConditionsImposition.ReducedVector(externalForces, elementsAssembly.BoundedDOFsVector);
            newSolver.LinearSolver = new LUFactorization();
            newSolver.ActivateNonLinearSolution = false;
            newSolver.SolveNewmark();
            //newSolver.PrintExplicitSolution();
            Results finalResults = new Results() { DynamicSolution = newSolver.explicitSolution, TimeSteps = newSolver.TimeAtEachStep, SelectedDOF = 1, SelectedInterval = 1, SolutionType = "Dynamic" };
            return finalResults;
        }
    }
}
