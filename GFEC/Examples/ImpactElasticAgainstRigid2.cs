using OpenTK.Graphics.ES11;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class ImpactElasticAgainstRigid2
    {
        //geometry & initial mesh
        private const double  width = 2.0;
        private const double height = 1.0;
        private const double thickness = 0.1;
        private const double xIntervals = 0.1;
        private const double yIntervals = 0.1;
        private const double gap = 0.045;

        //material
        private const double density = 800;
        private const double YoungMod = 200 * 1e9;

        private const int nodesInX = (int)(width/ xIntervals) + 1;
        private const int nodesInY = (int)(height / yIntervals) + 1;
        private const int nodesNumber = 2 * nodesInX * nodesInY;
        private const int elementsNumber = 2 * (nodesInX - 1) * (nodesInY - 1);
        private const int contacts = nodesInX;

        //external loads & boundary conditions
        private const double externalForce = 0;
        static int[] structuralBoundaryConditions;
        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;
        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int k = 1;
            for (int i = 0; i< nodesInY; i++)
            {
                for (int j = 0; j < nodesInX; j++)
                {
                    nodes[k] = new Node(j *xIntervals, i * yIntervals);// upper body nodes
                    k += 1;
                }
            }
            for (int i = 0; i < nodesInY; i++)
            {
                for (int j = 0; j < nodesInX; j++)
                {
                    nodes[k] = new Node(j * xIntervals, i * yIntervals - gap - height);// lower body nodes
                    k += 1;
                }
            }
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
            for (int i = 1; i <= nodesInX * (nodesInY - 1) + 1; i += nodesInX)
            {
                boundedDofs.Add(2 * i - 1); //left side support upper body
            }
            for (int i = nodesInX; i <= nodesInX * nodesInY; i += nodesInX)
            {
                boundedDofs.Add(2 * i - 1); //right side support upper body
            }
            for (int i = 1 + nodesNumber / 2; i <= nodesInX + nodesNumber / 2; i++)
            {
                if (i == 1 + nodesNumber / 2 ||
                    i == nodesInX + nodesNumber / 2)
                {
                    boundedDofs.Add(2 * i); //lower side support lower body
                }
                else
                {
                    boundedDofs.Add(2 * i - 1); //lower side support lower body
                    boundedDofs.Add(2 * i); //lower side support lower body
                }
            }
            for (int i = 1 + nodesNumber/2; i <= nodesInX * (nodesInY - 1) + 1 + nodesNumber / 2; i += nodesInX)
            {
                boundedDofs.Add(2 * i - 1); //left side support upper body
            }
            for (int i = nodesInX + nodesNumber / 2; i <= nodesNumber; i += nodesInX)
            {
                boundedDofs.Add(2 * i - 1); //right side support upper body
            }
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }
        private static double[] CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            for (int i = 1; i <= nodesInX; i++)
            {
                loadedStructuralDOFs.Add(i * 2);
            }
            externalForcesStructuralVector = new double[nodesNumber * 2];
            return externalForcesStructuralVector;
        }
        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int k = 1;
            for(int j = 0; j<nodesInY - 1; j++)
            {
                for (int i = 1; i < nodesInX; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + j * nodesInX }, { 2, i + 1 + j * nodesInX }, { 3, i + 1 + nodesInX + j * nodesInX }, { 4, i + nodesInX + j * nodesInX } };// upper body
                    k += 1;
                }
            }
            for (int j = 0; j < nodesInY - 1; j++)
            {
                for (int i = 1; i < nodesInX; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + j * nodesInX + nodesNumber/2 }, { 2, i + 1 + j * nodesInX + nodesNumber / 2 }, { 3, i + 1 + nodesInX + j * nodesInX + nodesNumber / 2 }, { 4, i + nodesInX + j * nodesInX + nodesNumber / 2 } };//lower body
                    k += 1;
                }
            }
            k = connectivity.Count + 1;
            for(int i = 1; i <= nodesInX; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, nodesNumber - nodesInX + i }, { 2, i} };//contacts
                k += 1;
            }
            return connectivity;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            double A = thickness * yIntervals;
            double contactElemArea = thickness * xIntervals;
            string type = "Quad4";
            string type2 = "ContactNtN2D";
            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= elementsNumber/2; i++)
            {
                elementProperties[i] = new ElementProperties(E, A, type);
            }
            for (int i = 1; i <= elementsNumber/2; i++)
            {
                elementProperties[i].Density = density;
                elementProperties[i].Thickness = thickness;
            }
            for (int i = elementsNumber / 2 + 1; i <= elementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(1000 * E, A, type);
            }
            for (int i = elementsNumber / 2 + 1; i <= elementsNumber; i++)
            {
                elementProperties[i].Density = density;
                elementProperties[i].Thickness = thickness;
            }
            for (int i = elementsNumber + 1; i <= elementsNumber + contacts; i++)
            {
                if (i == elementsNumber + 1 ||
                    i == elementsNumber + contacts)
                {
                    elementProperties[i] = new ElementProperties(E, contactElemArea / 2, type2);
                    elementProperties[i].Density = density;
                    elementProperties[i].Thickness = thickness;
                }
                else
                {
                    elementProperties[i] = new ElementProperties(E, contactElemArea, type2);
                    elementProperties[i].Density = density;
                    elementProperties[i].Thickness = thickness;
                }
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
            var AccelerationVector = new double[nodesNumber * 2];
            var DisplacementVector = new double[nodesNumber * 2];
            var VelocityVector = new double[nodesNumber * 2];
            //for (int i = 1; i < AccelerationVector.Length; i += 2)
            //{
            //    AccelerationVector[i] = -10;
            //}
            for (int i = 1; i < nodesNumber; i += 2)
            {
                VelocityVector[i] = -100.0;
            }
            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = BoundaryConditionsImposition.ReducedVector(AccelerationVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialDisplacementVector = BoundaryConditionsImposition.ReducedVector(DisplacementVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialVelocityVector = BoundaryConditionsImposition.ReducedVector(VelocityVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialTime = 0.0;
            //ExplicitSolver newSolver = new ExplicitSolver(0.0534, 100);
            ExplicitSolver newSolver = new ExplicitSolver(0.0005, 5000);
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
            newSolver.SolveExplicit();
            Tuple<Dictionary<int, double[]>, Dictionary<int, double>> solvectors = newSolver.GetResults();
            //int max = solvectors.Item1.OrderByDescending(m => m.Key).FirstOrDefault().Key;
            //elementsAssembly.UpdateDisplacements(solvectors.Item1.Single(m => m.Key == max).Value);
            Dictionary<int, double[]> allStepsSolutions = solvectors.Item1;
            //for (int i = 1; i < allStepsSolutions.Count; i++)
            for (int i = 0; i < allStepsSolutions.Count - 1; i++)
            {
                elementsAssembly.UpdateDisplacements(allStepsSolutions[i]);
            }
            ShowToGUI.PlotFinalGeometry(elementsAssembly);
            //double[] fullDynamicSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[93], elementsAssembly.BoundedDOFsVector);
            //double[] fullDynamicSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[94], elementsAssembly.BoundedDOFsVector);
            //double[] fullDynamicSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[95], elementsAssembly.BoundedDOFsVector);
            //double[] fullDynamicSol4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[96], elementsAssembly.BoundedDOFsVector);
            //double[] fullDynamicSol5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[97], elementsAssembly.BoundedDOFsVector);
            //VectorOperations.PrintVectorToFile(fullDynamicSol1, @"C:\Users\Public\Documents\Results1.dat");
            //VectorOperations.PrintVectorToFile(fullDynamicSol2, @"C:\Users\Public\Documents\Results2.dat");
            //VectorOperations.PrintVectorToFile(fullDynamicSol3, @"C:\Users\Public\Documents\Results3.dat");
            //VectorOperations.PrintVectorToFile(fullDynamicSol4, @"C:\Users\Public\Documents\Results4.dat");
            //VectorOperations.PrintVectorToFile(fullDynamicSol5, @"C:\Users\Public\Documents\Results5.dat");
            //newSolver.PrintExplicitSolution();
            Results finalResults = new Results() { DynamicSolution = newSolver.explicitSolution, TimeSteps = newSolver.TimeAtEachStep, SelectedDOF = 1, SelectedInterval = 1, SolutionType = "Dynamic" };
            return finalResults;
        }
    }
}
