using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public static class TwoThermalQuadsInContactExample
    {
        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            nodes[1] = new Node(0.0, 0.0);
            nodes[2] = new Node(1.0, 0.0);
            nodes[3] = new Node(1.0, 1.0);
            nodes[4] = new Node(0.0, 1.0);


            nodes[5] = new Node(0.0, 1.01);
            nodes[6] = new Node(1.0, 1.01);
            nodes[7] = new Node(1.0, 2.01);
            nodes[8] = new Node(0.0, 2.01);
            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {
            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            connectivity[1] = new Dictionary<int, int>() { { 1, 1 }, { 2, 2 }, { 3, 3 }, { 4, 4 } };
            connectivity[2] = new Dictionary<int, int>() { { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 8 } };
            connectivity[3] = new Dictionary<int, int>() { { 1, 4 }, { 2, 5 } };
            connectivity[4] = new Dictionary<int, int>() { { 1, 3 }, { 2, 6 } };

            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            nodeFAT[1] = new bool[] { true, false, false, false, false, false };
            nodeFAT[2] = new bool[] { true, false, false, false, false, false };
            nodeFAT[3] = new bool[] { true, false, false, false, false, false };
            nodeFAT[4] = new bool[] { true, false, false, false, false, false };
            nodeFAT[5] = new bool[] { true, false, false, false, false, false };
            nodeFAT[6] = new bool[] { true, false, false, false, false, false };
            nodeFAT[7] = new bool[] { true, false, false, false, false, false };
            nodeFAT[8] = new bool[] { true, false, false, false, false, false };
            return nodeFAT;
        }

        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double thermalCond = 1.0;
            double A = 1.0;
            string type = "Quad4Th";
            string type2 = "ContactNtN2DTh";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            elementProperties[1] = new ElementProperties();
            elementProperties[2] = new ElementProperties();
            elementProperties[3] = new ElementProperties();
            elementProperties[4] = new ElementProperties();
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
            assembly.BoundedDOFsVector = new int[] { 1, 2 };
            return assembly;
        }

        public static Results RunStaticExample()
        {
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();

            ISolver newSolu = new StaticSolver();
            newSolu.LinearScheme = new PCGSolver();
            newSolu.NonLinearScheme = new LoadControlledNewtonRaphson();
            newSolu.ActivateNonLinearSolver = true;
            newSolu.NonLinearScheme.numberOfLoadSteps = 10;

            double[] externalForces = new double[] { 0, 0, 0, 0, 250.0, 250.0 };
            newSolu.AssemblyData = elementsAssembly;
            newSolu.Solve(externalForces);
            newSolu.PrintSolution();
            double[] tempSolution = newSolu.GetSolution();
            return new Results();
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
