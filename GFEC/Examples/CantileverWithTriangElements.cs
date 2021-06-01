using OpenTK.Graphics.ES11;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class CantileverWithTriangElements
    {
        //geometry & initial mesh
        private const double length = 2.0;
        private const double height = 0.50;
        private const double thickness = 0.05;
        private const double xIntervals = 0.25/2;
        private const double yIntervals = 0.25/2;
        private const int  nodesInX = 21;
        private const int nodesInY = 5;


        //material
        private const double density = 8000;

        private const double YoungMod = 30.0 * 1e9;
        private const double PoissonRatio = 0.25;

        private const int nodesNumber = nodesInX * nodesInY;
        private const int elementsNumber = 2*(nodesInX - 1)* (nodesInY - 1);
        //external loads & boundary conditions
        private const double externalForce = -100000;
        static int[] structuralBoundaryConditions;
        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;
        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int k = 1;
            for(int i= 1; i <= nodesInX; i++)
            {
                nodes[k] = new Node((i-1)*xIntervals, 0);
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
            //for (int i = 1; i <= nodesInX; i++)
            //{
            //    nodes[k] = new Node((i - 1) * xIntervals, 5 * yIntervals);
            //    k += 1;
            //}
            //for (int i = 1; i <= nodesInX; i++)
            //{
            //    nodes[k] = new Node((i - 1) * xIntervals, 6 * yIntervals);
            //    k += 1;
            //}
            //for (int i = 1; i <= nodesInX; i++)
            //{
            //    nodes[k] = new Node((i - 1) * xIntervals, 7 * yIntervals);
            //    k += 1;
            //}
            //for (int i = 1; i <= nodesInX; i++)
            //{
            //    nodes[k] = new Node((i - 1) * xIntervals, 8 * yIntervals);
            //    k += 1;
            //}
            //for (int i = 1; i <= nodesInX; i++)
            //{
            //    nodes[k] = new Node((i - 1) * xIntervals, 9 * yIntervals);
            //    k += 1;
            //}
            //for (int i = 1; i <= nodesInX; i++)
            //{
            //    nodes[k] = new Node((i - 1) * xIntervals, 10 * yIntervals);
            //    k += 1;
            //}
            //for (int i = 1; i <= nodesInX; i++)
            //{
            //    nodes[k] = new Node((i - 1) * xIntervals, 11 * yIntervals);
            //    k += 1;
            //}
            //for (int i = 1; i <= nodesInX; i++)
            //{
            //    nodes[k] = new Node((i - 1) * xIntervals, 12 * yIntervals);
            //    k += 1;
            //}
            //for (int i = 1; i <= nodesInX; i++)
            //{
            //    nodes[k] = new Node((i - 1) * xIntervals, 13 * yIntervals);
            //    k += 1;
            //}
            //for (int i = 1; i <= nodesInX; i++)
            //{
            //    nodes[k] = new Node((i - 1) * xIntervals, 14 * yIntervals);
            //    k += 1;
            //}
            //for (int i = 1; i <= nodesInX; i++)
            //{
            //    nodes[k] = new Node((i - 1) * xIntervals, 15 * yIntervals);
            //    k += 1;
            //}
            //for (int i = 1; i <= nodesInX; i++)
            //{
            //    nodes[k] = new Node((i - 1) * xIntervals, 16 * yIntervals);
            //    k += 1;
            //}
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
            boundedDofs.Add(2 * 2*nodesInX - 1);
            boundedDofs.Add(2 * 2*nodesInX);
            boundedDofs.Add(2 * 3 * nodesInX - 1);
            boundedDofs.Add(2 * 3 * nodesInX);
            boundedDofs.Add(2 * 4 * nodesInX - 1);
            boundedDofs.Add(2 * 4 * nodesInX);
            boundedDofs.Add(2 * 5 * nodesInX - 1);
            boundedDofs.Add(2 * 5 * nodesInX);
            //boundedDofs.Add(2 * 6 * nodesInX - 1);
            //boundedDofs.Add(2 * 6 * nodesInX);
            //boundedDofs.Add(2 * 7 * nodesInX - 1);
            //boundedDofs.Add(2 * 7 * nodesInX);
            //boundedDofs.Add(2 * 8 * nodesInX - 1);
            //boundedDofs.Add(2 * 8 * nodesInX);
            //boundedDofs.Add(2 * 9 * nodesInX - 1);
            //boundedDofs.Add(2 * 9 * nodesInX);
            //boundedDofs.Add(2 *10 *  nodesInX - 1);
            //boundedDofs.Add(2 *10 * nodesInX);
            //boundedDofs.Add(2 * 11 * nodesInX - 1);
            //boundedDofs.Add(2 * 11 * nodesInX);
            //boundedDofs.Add(2 * 12 * nodesInX - 1);
            //boundedDofs.Add(2 * 12 * nodesInX);
            //boundedDofs.Add(2 * 13 * nodesInX - 1);
            //boundedDofs.Add(2 * 13 * nodesInX);
            //boundedDofs.Add(2 * 14 * nodesInX - 1);
            //boundedDofs.Add(2 * 14 * nodesInX);
            //boundedDofs.Add(2 * 15 * nodesInX - 1);
            //boundedDofs.Add(2 * 15 * nodesInX);
            //boundedDofs.Add(2 * 16 * nodesInX - 1);
            //boundedDofs.Add(2 * 16 * nodesInX);
            //boundedDofs.Add(2 * 17 * nodesInX - 1);
            //boundedDofs.Add(2 * 17 * nodesInX);
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }
        private static double[] CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            //loadedStructuralDOFs.Add(1);
            //loadedStructuralDOFs.Add((1 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((2 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((3 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((4 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((5 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((6 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((7 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((8 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((9 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((10 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((11 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((12 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((13 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((14 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((15 * nodesInX + 1) * 2 - 1);
            //loadedStructuralDOFs.Add((16 * nodesInX + 1) * 2 - 1);
            loadedStructuralDOFs.Add((4 * nodesInX + 1) * 2);

            //loadedStructuralDOFs.Add(1);
            //loadedStructuralDOFs.Add(18 * 2 - 1);
            //loadedStructuralDOFs.Add(35 * 2 - 1);
            //loadedStructuralDOFs.Add(52 * 2 - 1);
            //loadedStructuralDOFs.Add(69 * 2 - 1);
            externalForcesStructuralVector = new double[nodesNumber * 2];
            //externalForcesStructuralVector = new double[nodesInY * nodesInY];
            return externalForcesStructuralVector;
        }
        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int k = 1;
            for (int i = 1; i < nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            for (int i = 1 + nodesInX; i < 2 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            for (int i = 1 + 2 * nodesInX; i < 3 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            for (int i = 1 + 3 * nodesInX; i < 4 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            for (int i = 1 + 4 * nodesInX; i < 5 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            for (int i = 1 + 5 * nodesInX; i < 6 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            for (int i = 1 + 6 * nodesInX; i < 7 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            for (int i = 1 + 7 * nodesInX; i < 8 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }

            for (int i = 1 + 8 * nodesInX; i < 9 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            for (int i = 1 + 9 * nodesInX; i < 10 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            for (int i = 1 + 10 * nodesInX; i < 11 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            for (int i = 1 + 11 * nodesInX; i < 12 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            for (int i = 1 + 12 * nodesInX; i < 13 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            for (int i = 1 + 13 * nodesInX; i < 14 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            for (int i = 1 + 14 * nodesInX; i < 15 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            for (int i = 1 + 15 * nodesInX; i < 16 * nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + nodesInX } };
                k += 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 1 }, { 2, i + 1 + nodesInX }, { 3, i + nodesInX } };
                k += 1;
            }
            //int k = 1;
            //for (int i = 1; i < nodesInX; i++)
            //{
            //    connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + 1 + nodesInX }, { 4, i + nodesInX } };
            //    k += 1;
            //}
            //for (int i = 1 + nodesInX; i < 2 * nodesInX; i++)
            //{
            //    connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + 1 + nodesInX }, { 4, i + nodesInX } };
            //    k += 1;
            //}
            //for (int i = 1 + 2 * nodesInX; i < 3 * nodesInX; i++)
            //{
            //    connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + 1 + nodesInX }, { 4, i + nodesInX } };
            //    k += 1;
            //}
            //for (int i = 1 + 3 * nodesInX; i < 4 * nodesInX; i++)
            //{
            //    connectivity[k] = new Dictionary<int, int>() { { 1, i }, { 2, i + 1 }, { 3, i + 1 + nodesInX }, { 4, i + nodesInX } };
            //    k += 1;
            //}
            return connectivity;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            double A = thickness * yIntervals;
            string type = "Triangle3";
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
            //newSolver.NonLinearScheme = new LoadControlledNewtonRaphson();
            //newSolver.NonLinearScheme.numberOfLoadSteps = 10;
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
