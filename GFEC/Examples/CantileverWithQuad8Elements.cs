using OpenTK.Graphics.ES11;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class CantileverWithQuad8Elements
    {
        //geometry & initial mesh
        private const double length = 2.0;
        private const double height = 0.50;
        private const double thickness = 0.05;
        private const double xIntervals = 0.25 / 2;
        private const double yIntervals = 0.25 / 2;
        private const int nodesInX = 17;
        private const int nodesInY = 5;
        private const int addedNodes = (2 * nodesInX - 1) *(nodesInY - 1) + nodesInX - 1;


        //material
        private const double density = 8000;

        private const double YoungMod = 30.0 * 1e9;
        private const double PoissonRatio = 0.25;

        private const int nodesNumber = nodesInX * nodesInY + addedNodes;
        private const int elementsNumber = 2 * (nodesInX - 1) * (nodesInY - 1);
        //external loads & boundary conditions
        private const double externalForce = -100000;
        static int[] structuralBoundaryConditions;
        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;
        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int k = 1;
            for (int i = 1; i <= nodesInX; i++)
            {
                nodes[k] = new Node((i - 1) * xIntervals, 0);
                k += 1;
            }
            for (int i = 1; i <= nodesInX; i++)
            {
                nodes[k] = new Node((i - 1) * xIntervals, yIntervals);
                k += 1;
            }
            for (int i = 1; i <= nodesInX; i++)
            {
                nodes[k] = new Node((i - 1) * xIntervals, 2 * yIntervals);
                k += 1;
            }
            for (int i = 1; i <= nodesInX; i++)
            {
                nodes[k] = new Node((i - 1) * xIntervals, 3 * yIntervals);
                k += 1;
            }
            for (int i = 1; i <= nodesInX; i++)
            {
                nodes[k] = new Node((i - 1) * xIntervals, 4 * yIntervals);
                k += 1;
            }
            for(int j = 1; j<= 2 * nodesInY - 1; j++)
            {
                for(int i = 1; i<nodesInX; i++)
                {
                    if(j%2 != 0)
                    {
                        nodes[k] = new Node((i - 1) * xIntervals + xIntervals/2.0, (j - 1) * yIntervals/2);
                        k += 1;
                    }
                    else
                    {
                        nodes[k] = new Node((i - 1) * xIntervals, (j - 1) * yIntervals / 2);
                        k += 1;
                        if (i == nodesInX - 1)
                        {
                            nodes[k] = new Node((nodesInX - 1) * xIntervals, (j - 1) * yIntervals / 2);
                            k += 1;
                        }
                    }
                }
            }
            return nodes;
        }
        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            //for (int i = 1; i <= nodesInY * nodesInY; i++)
            for (int i = 1; i <= nodesNumber; i++)
            {
                nodeFAT[i] = new bool[] { true, true, false, false, false, false };
            }
            return nodeFAT;
        }
        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            boundedDofs.Add(2 * nodesInX - 1);
            boundedDofs.Add(2 * nodesInX);
            boundedDofs.Add(2 * 2 * nodesInX - 1);
            boundedDofs.Add(2 * 2 * nodesInX);
            boundedDofs.Add(2 * 3 * nodesInX - 1);
            boundedDofs.Add(2 * 3 * nodesInX);
            boundedDofs.Add(2 * 4 * nodesInX - 1);
            boundedDofs.Add(2 * 4 * nodesInX);
            boundedDofs.Add(2 * 5 * nodesInX - 1);
            boundedDofs.Add(2 * 5 * nodesInX);
            boundedDofs.Add(2 * (nodesInY * nodesInX + nodesInX - 1) - 1);
            boundedDofs.Add(2 * (nodesInY * nodesInX + nodesInX - 1));
            boundedDofs.Add(2 * (nodesInY * nodesInX + 2 * nodesInX - 1) - 1);
            boundedDofs.Add(2 * (nodesInY * nodesInX + 2 * nodesInX - 1));
            boundedDofs.Add(2 * (nodesInY * nodesInX + 3 * nodesInX - 2) - 1);
            boundedDofs.Add(2 * (nodesInY * nodesInX + 3 * nodesInX - 2));
            boundedDofs.Add(2 * (nodesInY * nodesInX + 4 * nodesInX - 2) - 1);
            boundedDofs.Add(2 * (nodesInY * nodesInX + 4 * nodesInX - 2));
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }
        private static double[] CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            loadedStructuralDOFs.Add((2 * nodesInX + 1) * 2);
            externalForcesStructuralVector = new double[nodesNumber * 2];
            return externalForcesStructuralVector;
        }
        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int k = 1;
            for (int i = 1; i < nodesInX; i++)
            {
                int second = nodesInX * nodesInY + i;
                int fourth = nodesInX * nodesInY + (nodesInX - 1) + i + 1;
                int sixth = nodesInX * nodesInY + (2*nodesInX - 1) + i;
                int last = fourth - 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, second }, { 3, i + 1},
                    {4, fourth}, {5, i + nodesInX + 1}, {6, sixth}, {7, i + nodesInX}, {8, last}};
                k += 1;
            }
            for (int i = 1 + nodesInX; i < 2 * nodesInX; i++)
            {
                int second = nodesInX * nodesInY + i + (nodesInX - 1);
                int fourth = nodesInX * nodesInY + (nodesInX - 1) + i + 1 + (nodesInX - 1);
                int sixth = nodesInX * nodesInY + (2 * nodesInX - 1) + i + (nodesInX - 1);
                int last = fourth - 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, second }, { 3, i + 1},
                    {4, fourth}, {5, i + nodesInX + 1}, {6, sixth}, {7, i + nodesInX}, {8, last}};
                k += 1;
            }
            for (int i = 1 + 2 * nodesInX; i < 3 * nodesInX; i++)
            {
                int second = nodesInX * nodesInY + i + 2*(nodesInX - 1);
                int fourth = nodesInX * nodesInY + 2*(nodesInX - 1) + i + 1 + (nodesInX - 1);
                int sixth = nodesInX * nodesInY + (2 * nodesInX - 1) + i + 2*(nodesInX - 1);
                int last = fourth - 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, second }, { 3, i + 1},
                    {4, fourth}, {5, i + nodesInX + 1}, {6, sixth}, {7, i + nodesInX}, {8, last}};
                k += 1;
            }
            for (int i = 1 + 3 * nodesInX; i < 4 * nodesInX; i++)
            {
                int second = nodesInX * nodesInY + i + 3 * (nodesInX - 1);
                int fourth = nodesInX * nodesInY + 3 * (nodesInX - 1) + i + 1 + (nodesInX - 1);
                int sixth = nodesInX * nodesInY + (2 * nodesInX - 1) + i + 3 * (nodesInX - 1);
                int last = fourth - 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, second }, { 3, i + 1},
                    {4, fourth}, {5, i + nodesInX + 1}, {6, sixth}, {7, i + nodesInX}, {8, last}};
                k += 1;
            }
            return connectivity;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            double A = thickness * yIntervals;
            string type = "Quad8";
            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= elementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E, PoissonRatio, A, thickness, density, type);
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
            ShowToGUI.PlotInitialGeometry(elementsAssembly);
            ISolver newSolver = new StaticSolver();
            newSolver.LinearScheme = new LUFactorization();
            newSolver.ActivateNonLinearSolver = false;
            newSolver.AssemblyData = elementsAssembly;
            double[] externalForces = externalForcesStructuralVector;
            foreach (var dof in loadedStructuralDOFs)
            {
                externalForces[dof - 1] = externalForce;
                //if (dof == 1 || dof == (16 * nodesInX + 1) * 2 - 1)
                //{
                //    externalForces[dof - 1] = externalForce / 2;
                //}
                //else
                //{
                //    externalForces[dof - 1] = externalForce;
                //}
            }
            newSolver.Solve(BoundaryConditionsImposition.ReducedVector(externalForces, elementsAssembly.BoundedDOFsVector));
            double[] solvector = newSolver.GetSolution();
            elementsAssembly.UpdateDisplacements(solvector);
            ShowToGUI.PlotFinalGeometry(elementsAssembly);
            List<double[]> solutions = new List<double[]>();
            //newSolver.PrintExplicitSolution();
            return new Results() { NonlinearSolution = solutions, SelectedDOF = 2, SolutionType = "Nonlinear" };

        }
    }
}
