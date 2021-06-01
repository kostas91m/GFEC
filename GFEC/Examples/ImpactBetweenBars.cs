using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class ImpactBetweenBars
    {
        private const double length = 5.0;
        private const double gap = 0.00000001;
        private const double sectionArea = 0.315;
        private const double density = 800;
        private const int nodesNumber = 4;
        private const double YoungMod = 200 * 1e9;
        private const int elementsNumber = 2;
        private const double externalForce = 50 * 1e7;
        
        static int[] structuralBoundaryConditions;
        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;

        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            nodes[1] = new Node(0, 0);
            nodes[2] = new Node(length, 0);
            nodes[3] = new Node(length + gap, 0);
            nodes[4] = new Node(2 * length + gap, 0);
            return nodes;
        }
        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= nodesNumber; i++)
            {
                nodeFAT[i] = new bool[] { true, true, false, false, false, false };
            }
            return nodeFAT;
        }
        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            boundedDofs.Add(1); //right side support 
            boundedDofs.Add(2); //right side support 
            boundedDofs.Add(7); //left side support 
            boundedDofs.Add(8); //left side support 
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }
        private static double[] CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            loadedStructuralDOFs.Add(3);
            externalForcesStructuralVector = new double[nodesNumber * 2];
            return externalForcesStructuralVector;
        }
        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            connectivity[1] = new Dictionary<int, int>() { { 1, 1 }, { 2, 2 } };
            connectivity[2] = new Dictionary<int, int>() { { 1, 3 }, { 2, 4 } };
            connectivity[3] = new Dictionary<int, int>() { { 1, 2 }, { 2, 3 } };//contact
            return connectivity;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            double A = sectionArea;
            string type = "Bar2D";
            string type2 = "ContactNtN2D";
            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= elementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E, A, type);
                elementProperties[i].Density = density;
            }
            elementProperties[elementsNumber + 1] = new ElementProperties(E, A, type2);
            elementProperties[elementsNumber + 1].Density = density;
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
            ShowToGUI.PlotInitialGeometry(elementsAssembly);
            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = new double[] { 10.0, 0.0, 10.0, 0.0 };
            initialValues.InitialDisplacementVector = new double[] { 0.0, 0.0, 0.0, 0.0 };
            initialValues.InitialVelocityVector = new double[] { 0.0, 0.0, 0.0, 0.0 };
            initialValues.InitialTime = 0.0;

            ExplicitSolver newSolver = new ExplicitSolver(1, 100);
            newSolver.Assembler = elementsAssembly;
            double[] externalForces = externalForcesStructuralVector;
            foreach (var dof in loadedStructuralDOFs)
            {
                externalForces[dof - 1] = externalForce;

            }
            newSolver.InitialValues = initialValues;
            newSolver.ExternalForcesVector = BoundaryConditionsImposition.ReducedVector(externalForces, elementsAssembly.BoundedDOFsVector);
            newSolver.LinearSolver = new LUFactorization();
            newSolver.ActivateNonLinearSolution = true;
            newSolver.SolveNewmark();
            Tuple<Dictionary<int, double[]>, Dictionary<int, double>> solvectors = newSolver.GetResults();
            int max = solvectors.Item1.OrderByDescending(m => m.Key).FirstOrDefault().Key;
            elementsAssembly.UpdateDisplacements(solvectors.Item1.Single(m => m.Key == max).Value);
            ShowToGUI.PlotFinalGeometry(elementsAssembly);
            //newSolver.PrintExplicitSolution();
            Results finalResults = new Results() { DynamicSolution = newSolver.explicitSolution, TimeSteps = newSolver.TimeAtEachStep, SelectedDOF = 1, SelectedInterval = 1, SolutionType = "Dynamic" };
            return finalResults;
        }
    }
}
