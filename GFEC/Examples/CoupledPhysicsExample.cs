using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public static class CoupledPhysicsExample
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
            nodeFAT[1] = new bool[] { true, true, false, false, false, false };
            nodeFAT[2] = new bool[] { true, true, false, false, false, false };
            nodeFAT[3] = new bool[] { true, true, false, false, false, false };
            nodeFAT[4] = new bool[] { true, true, false, false, false, false };
            nodeFAT[5] = new bool[] { true, true, false, false, false, false };
            nodeFAT[6] = new bool[] { true, true, false, false, false, false };
            nodeFAT[7] = new bool[] { true, true, false, false, false, false };
            nodeFAT[8] = new bool[] { true, true, false, false, false, false };
            return nodeFAT;
        }

        private static Dictionary<int, bool[]> CreateThermalNodeFAT()
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
            double E = 200.0e9;
            double A = 0.01;
            string type = "Quad4";
            string type2 = "ContactNtN2D";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            elementProperties[1] = new ElementProperties(E, A, type);
            elementProperties[2] = new ElementProperties(E, A, type);
            elementProperties[3] = new ElementProperties(E / 1000.0, A, type2);
            elementProperties[4] = new ElementProperties(E / 1000.0, A, type2);
            for (int i = 1; i <= 2; i++)
            {
                elementProperties[i].Density = 8000.0;
                elementProperties[i].Thickness = 0.1;
            }

            return elementProperties;
        }

        private static Dictionary<int, IElementProperties> CreateThermalElementProperties()
        {
            double thermalCond = 60.5;
            double A = 0.5;
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
            assembly.BoundedDOFsVector = new int[] { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15 };
            return assembly;
        }

        private static IAssembly CreateThermalAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes();
            assembly.ElementsConnectivity = CreateConnectivity();
            assembly.ElementsProperties = CreateThermalElementProperties();
            assembly.NodeFreedomAllocationList = CreateThermalNodeFAT();
            assembly.BoundedDOFsVector = new int[] { 1, 2 };
            return assembly;
        }

        public static Results RunStaticExample()
        {
            ///Structural settings
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();

            ISolver newSolu = new StaticSolver();
            newSolu.LinearScheme = new PCGSolver();
            newSolu.NonLinearScheme = new LoadControlledNewtonRaphson();
            newSolu.ActivateNonLinearSolver = true;
            newSolu.NonLinearScheme.numberOfLoadSteps = 1;


            //Thermal Settings
            IAssembly elementsAssembly2 = CreateThermalAssembly();
            elementsAssembly2.CreateElementsAssembly();
            elementsAssembly2.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix2 = elementsAssembly2.CreateTotalStiffnessMatrix();

            ISolver thermalSolution = new StaticSolver();
            thermalSolution.LinearScheme = new CholeskyFactorization();
            thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
            thermalSolution.ActivateNonLinearSolver = true;
            thermalSolution.NonLinearScheme.numberOfLoadSteps = 10;

            //double[] externalFlux = new double[] { 0, 0, 0, 0, 250.0, 250.0 };
            //newSolu2.AssemblyData = elementsAssembly;
            //newSolu2.Solve(externalForces);
            //newSolu2.PrintSolution();
            //double[] tempSolution = newSolu.GetSolution();
        
            double[] solVector2 = new double[6];
            List<double[]> structuralSolutions = new List<double[]>();
            for (int i = 1; i <= 5; i++)
            {
                newSolu.NonLinearScheme = new LoadControlledNewtonRaphson(solVector2);
                double[] externalForces2 = new double[] { 0, 0, 0, 0, -10000.0*i, -10000.0*i };
                newSolu.AssemblyData = elementsAssembly;
                newSolu.Solve(externalForces2);
                solVector2 = newSolu.GetSolution();
                structuralSolutions.Add(solVector2);
            }
            Dictionary<int, double[]> intForces = newSolu.GetInternalForces();
            Dictionary<int, double[]> elementInternalForces = elementsAssembly.GetElementsInternalForces(structuralSolutions[0]);
            List<string> elementTypes = elementsAssembly.GetElementsType();

            double[] temperatures = new double[6];
            List<double[]> thermalSolutions = new List<double[]>();
            for (int i = 1; i <= 5; i++)
            {
                thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson(temperatures);
                double[] externalHeatFlux = new double[] { 0, 0, 0, 0, 250.0, 250.0 };
                thermalSolution.AssemblyData = elementsAssembly2;
                thermalSolution.Solve(externalHeatFlux);
                temperatures = thermalSolution.GetSolution();
                thermalSolutions.Add(temperatures);
            }


            //double[] completeFinalSolutionVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solVector2, new int[] { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15 });            Dictionary<int, INode> finalNodesList = new Dictionary<int, INode>();
            //finalNodesList = Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, completeFinalSolutionVector);
             
            return new Results() { NonlinearSolution = structuralSolutions, SelectedDOF=2, SolutionType="Nonlinear" };
        }

        private static void NonLinearScheme_convergenceResult(object sender, string e)
        {
            throw new NotImplementedException();
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
