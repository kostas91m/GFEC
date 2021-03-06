﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public static class LinearTrussExample
    {
        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            nodes[1] = new Node(0.0, 0.0);
            nodes[2] = new Node(40.0, 0.0);
            nodes[3] = new Node(40.0, 30.0);
            nodes[4] = new Node(0.0, 30.0);
            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {
            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            connectivity[1] = new Dictionary<int, int>() { { 1, 1 }, { 2, 2 } };
            connectivity[2] = new Dictionary<int, int>() { { 1, 3 }, { 2, 2 } };
            connectivity[3] = new Dictionary<int, int>() { { 1, 1 }, { 2, 3 } };
            connectivity[4] = new Dictionary<int, int>() { { 1, 4 }, { 2, 3 } };
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            nodeFAT[1] = new bool[] { true, true, false, false, false, false };
            nodeFAT[2] = new bool[] { true, true, false, false, false, false };
            nodeFAT[3] = new bool[] { true, true, false, false, false, false };
            nodeFAT[4] = new bool[] { true, true, false, false, false, false };
            return nodeFAT;
        }

        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = 29.5e6;
            double A = 1;
            string type = "Bar2D";
            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            elementProperties[1] = new ElementProperties(E, A, type);
            elementProperties[2] = new ElementProperties(E, A, type);
            elementProperties[3] = new ElementProperties(E, A, type);
            elementProperties[4] = new ElementProperties(E, A, type);
            return elementProperties;
        }

        private static IAssembly CreateAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes();
            assembly.ElementsConnectivity = CreateConnectivity();
            assembly.ElementsProperties = CreateElementProperties();
            assembly.NodeFreedomAllocationList = CreateNodeFAT();
            assembly.BoundedDOFsVector = new int[] { 1, 2, 4, 7, 8 };
            return assembly;
        }

        public static Results RunExample()
        {
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();

            ISolver newSolu = new StaticSolver();
            newSolu.LinearScheme = new PCGSolver();

            double[] externalForces = new double[] { 20000, 0, -25000 };
            newSolu.AssemblyData = elementsAssembly;
            newSolu.Solve(externalForces);
            newSolu.PrintSolution();
            return new Results();
        }
    }
}

