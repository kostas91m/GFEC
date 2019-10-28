using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public static class ThermalExample
    {
        private static Dictionary<int, INode> CreateNodes()
        {
            double a = 1.0;

            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            nodes[1] = new Node(0.0, 0.0);
            nodes[2] = new Node(a, 0.0);
            nodes[3] = new Node(2.0 * a, 0.0);
            nodes[4] = new Node(3.0 * a, 0.0);

            nodes[5] = new Node(0.0, a);
            nodes[6] = new Node(a, a);
            nodes[7] = new Node(2.0 * a, a);
            nodes[8] = new Node(3.0 * a, a);

            nodes[9] = new Node(0.0, 2.0 * a);
            nodes[10] = new Node(a, 2.0 * a);
            nodes[11] = new Node(2.0 * a, 2.0 * a);
            nodes[12] = new Node(3.0 * a, 2.0 * a);

            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {
            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            connectivity[1] = new Dictionary<int, int>() { { 1, 1 }, { 2, 2 }, { 3, 6 }, { 4, 5 } };
            connectivity[2] = new Dictionary<int, int>() { { 1, 2 }, { 2, 3 }, { 3, 7 }, { 4, 6 } };
            connectivity[3] = new Dictionary<int, int>() { { 1, 3 }, { 2, 4 }, { 3, 8 }, { 4, 7 } };
            connectivity[4] = new Dictionary<int, int>() { { 1, 5 }, { 2, 6 }, { 3, 10 }, { 4, 9 } };
            connectivity[5] = new Dictionary<int, int>() { { 1, 6 }, { 2, 7 }, { 3, 11 }, { 4, 10 } };
            connectivity[6] = new Dictionary<int, int>() { { 1, 7 }, { 2, 8 }, { 3, 12 }, { 4, 11 } };

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
            nodeFAT[9] = new bool[] { true, false, false, false, false, false };
            nodeFAT[10] = new bool[] { true, false, false, false, false, false };
            nodeFAT[11] = new bool[] { true, false, false, false, false, false };
            nodeFAT[12] = new bool[] { true, false, false, false, false, false };
            return nodeFAT;
        }

        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double thermaConductivity = 1.0;
            string type = "Quad4Th";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= 6; i++)
            {
                elementProperties[i] = new ElementProperties();
            }
            foreach (var element in elementProperties)
            {
                element.Value.ElementType = type;
                element.Value.ThermalConductivity = thermaConductivity;
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
            assembly.BoundedDOFsVector = new int[] { 4, 8, 9, 10, 11, 12 };
            return assembly;
        }

        public static Results RunStaticExample()
        {
            double T0 = 1.0;
            double kc = 1.0;

            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();

            ISolver newSolu = new StaticSolver();
            newSolu.LinearScheme = new PCGSolver();
            newSolu.ActivateNonLinearSolver = false;

            double[] externalForces = new double[] { 0.0, 0.0, 0.0, (T0 + Math.Sqrt(3.0) * T0) * kc / 6.0, (2.0 * T0 + Math.Sqrt(3.0) * T0 + T0) * kc / 6.0, (T0 + Math.Sqrt(3.0) * T0) * kc / 6.0 };
            newSolu.AssemblyData = elementsAssembly;
            newSolu.Solve(externalForces);
            newSolu.PrintSolution();
            double[] kati = newSolu.GetSolution();
            return new Results();
        }
    }
}
