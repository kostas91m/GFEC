using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public class  NewExampleContacts
    {
        public static ISolver structuralSolution;
        private const double thickness = 0.1;
        private const int nodesCount = 2 * 126;
        private const int nodesInX = 21;
        private const int addedNodes = 2;
        private const int nodesInY = 6;
        private const double xInterv = 0.1;
        private const double yInterv = 0.1;
        private const double gap = 0.01;
        private const int elementsCount = 2 * 100;
        private const int contactElementsCount = nodesInX;
        private const double YoungMod = 200 * 1e9;
        static double[] nodalForcesVector;
        static int[] boundaryConditions;
        static List<int> LoadedStructuralDOFs = new List<int>();
        private static void CreateBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            for(int i = 1; i <= nodesInY; i++)
            {
                //boundedDofs.Add(((i - 1) * nodesInX + i) * 2);
                boundedDofs.Add(((i - 1) * nodesInX + i) * 2 - 1);
            }
            for (int i = 1; i <= nodesInY; i++)
            {
                //boundedDofs.Add(((i - 1) * nodesInX + i)*2);
                boundedDofs.Add((i * nodesInX) * 2 - 1);
            }
            for (int i = 2; i <= nodesInY; i++)
            {
                //boundedDofs.Add(((i - 1) * nodesInX + i)*2);
                boundedDofs.Add(((i - 1) * nodesInX + i) * 2 - 1 + nodesCount);
            }
            for (int i = 2; i <= nodesInY; i++)
            {
                //boundedDofs.Add(((i - 1) * nodesInX + i)*2);
                boundedDofs.Add((i * nodesInX) * 2 - 1 + nodesCount);
            }
            for (int i = 1; i <= nodesInX; i++)
            {
                boundedDofs.Add(i*2 + nodesCount);
                boundedDofs.Add(i*2 + nodesCount - 1);
            }
            for(int i = nodesCount + 1; i<=nodesCount + addedNodes; i++)
            {
                boundedDofs.Add(i* 2);
                boundedDofs.Add(i * 2 - 1);
            }
            boundaryConditions = boundedDofs.ToArray<int>();
        }
        private static void CreateLoadVector()
        {
            List<int> loadedDOFs = new List<int>(); 
            for(int i = nodesCount/2 - 1; i>=nodesCount/2 - nodesInX + 2; i-=1)
            {
                loadedDOFs.Add(2 * i);
            }
            LoadedStructuralDOFs = loadedDOFs;
            nodalForcesVector = new double[(nodesCount + addedNodes) * 2];
        }
        private static Dictionary <int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int k = 1;
            for(int i = 1; i<= nodesInY; i++)
            {
                for(int j = 1; j<=nodesInX; j++)
                {
                    Node node = new Node((j - 1) * xInterv, (i - 1) * yInterv);
                    nodes.Add(k, node);
                    k += 1;
                }
            }
            for (int i = 1; i <= nodesInY; i++)
            {
                for (int j = 1; j <= nodesInX; j++)
                {
                    Node node = new Node((j - 1) * xInterv, (i - 1) * yInterv - gap - (nodesInY - 1) * (yInterv));
                    nodes.Add(k, node);
                    k += 1;
                }
            }
            Node node2 = new Node(0, nodesInY * (yInterv));
            nodes.Add(k, node2);
            k += 1;
            Node node3 = new Node((nodesInX - 1) * xInterv, nodesInY * (yInterv));
            nodes.Add(k, node3);
            return nodes;
        }
        private static Dictionary<int, Dictionary<int,int>> CreateConnectivity()
        {
            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int k = 1;
            for (int j = 0; j <= nodesInY - 2; j++)
            {
                for (int i = 1; i <= nodesInX - 1; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + j * nodesInX }, { 2, i + 1 + j * nodesInX }, { 3, i + 1 + nodesInX + j * nodesInX }, { 4, i + nodesInX + j * nodesInX } };
                    k += 1;
                }

            }
            for (int j = 0; j <= nodesInY - 2; j++)
            {
                for (int i = 1; i <= contactElementsCount - 1; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, nodesInX * nodesInY +  i + j * nodesInX }, { 2, nodesInX * nodesInY + i + 1 + j * nodesInX }, { 3, nodesInX * nodesInY + i + 1 + nodesInX + j * nodesInX }, { 4, nodesInX * nodesInY + i + nodesInX + j * nodesInX } };
                    k += 1;
                }

            }
            //contact elements
            k = connectivity.Count + 1;
            for(int i = 1; i<=nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, nodesCount - nodesInX + i }};
                k += 1;
            }
            connectivity[k] = new Dictionary<int, int>() { { 1, nodesInX * nodesInY - nodesInX + 1}, { 2, nodesCount + 1 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, nodesInX * nodesInY }, { 2, nodesCount + 2} };
            return connectivity;
        }
        private static Dictionary <int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for(int i = 1; i<= nodesCount + addedNodes; i++)
            {
                nodeFAT[i] = new bool[] { true, true, false, false, false, false};
            }
            return nodeFAT;
        }
        private static Dictionary<int, IElementProperties > CreateElementProperties()
        {
            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            string type = "Quad4";
            string type2 = "ContactNtN2D";
            string type3 = "Bar2D";
            double density = new double();
            double A = xInterv * yInterv;
            double contactArea = xInterv * thickness;


            for (int i = 1; i<= elementsCount; i++)
            {
                elementProperties[i] = new ElementProperties(YoungMod, A, type);
            }
            for (int i = 1; i <= elementsCount; i++)
            {
                elementProperties[i].Density = density;
                elementProperties[i].Thickness = thickness;

            }
            for (int i = elementsCount + 1; i <= elementsCount + contactElementsCount; i++)
            {
                if (i == elementsCount + 1 || i == elementsCount + contactElementsCount)
                {
                    elementProperties[i] = new ElementProperties(YoungMod, contactArea/2, type2);
                    elementProperties[i].Thickness = thickness;
                }
                else
                {
                    elementProperties[i] = new ElementProperties(YoungMod, contactArea, type2);
                    elementProperties[i].Thickness = thickness;
                }
            }
            int count = elementProperties.Count;
            elementProperties[count + 1] = new ElementProperties(YoungMod/100, A/2, type3);
            elementProperties[count + 2] = new ElementProperties(YoungMod/100, A/2, type3);
            return elementProperties;
        }
        private static IAssembly CreateAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes();
            assembly.ElementsConnectivity = CreateConnectivity();
            assembly.ElementsProperties = CreateElementProperties();
            assembly.NodeFreedomAllocationList = CreateNodeFAT();
            CreateBoundaryConditions();
            CreateLoadVector();
            assembly.BoundedDOFsVector = boundaryConditions;
            return assembly;
        }

        public static Results RunStaticExample()
        {
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();
            ShowToGUI.PlotInitialGeometry(elementsAssembly);
            double[] externalForces = nodalForcesVector;
            foreach(var dof in LoadedStructuralDOFs)
            {
                externalForces[dof - 1] = -4 * 100 * 1e3;
            }
            double[] reducedNodalForces = BoundaryConditionsImposition.ReducedVector(externalForces, elementsAssembly.BoundedDOFsVector);
            structuralSolution.AssemblyData = elementsAssembly;
            structuralSolution.LinearScheme = new LUFactorization();
            structuralSolution.ActivateNonLinearSolver = true;
            //structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 40;
            structuralSolution.NonLinearScheme.Tolerance = 1e-4;
            structuralSolution.Solve(reducedNodalForces);
            var solution = structuralSolution.GetSolution();
            elementsAssembly.UpdateDisplacements(solution);
            ShowToGUI.PlotFinalGeometry(elementsAssembly);
            Dictionary<int, double[]> allstepssolutions = structuralSolution.GetAllStepsSolutions();
            List<double[]> solutions = new List<double[]>();
            structuralSolution.PrintSolution();
            return new Results() { NonlinearSolution = solutions, SelectedDOF = 2, SolutionType = "Nonlinear" };
        }

    }
}
