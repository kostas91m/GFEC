using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public static class ThreeTrusses
    {

            public static ISolver structuralSolution;
            static int[] structuralBoundaryConditions;
            const double length = 5.0;
            const double gap = 0.1;

        //External loads
        const double externalStructuralLoad = 50.0;

            static List<int> loadedStructuralDOFs; 
            static double[] externalForcesStructuralVector; 

            const double YoungMod = 1.0 * 1e3;
            const double density = 8000.0;
            const double area = 1.0;



            private static void CreateStructuralBoundaryConditions()
            {
            List<int> boundedDofs = new List<int>();
            boundedDofs.Add(1);
            boundedDofs.Add(2);
            boundedDofs.Add(4);
            boundedDofs.Add(6);
            boundedDofs.Add(8);
            boundedDofs.Add(9);
            boundedDofs.Add(10);
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
            }

            private static void CreateStructuralLoadVector()
            {
                loadedStructuralDOFs = new List<int>();
                loadedStructuralDOFs.Add(3);
                externalForcesStructuralVector = new double[5 * 2];
            }

            private static Dictionary<int, INode> CreateNodes()
            {
            
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int k;
            k = 1;
            for (int j = 0; j <= 2; j++)
            {
            nodes[k] = new Node(j * length, 0);
            k += 1;
            }
            nodes[k] = new Node(2 * length + gap, 0);
            k += 1;
            nodes[k] = new Node(3 * length + gap, 0);
            return nodes;
            }

            private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
            {

                Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
                int k = 1;
                for (int j = 1; j < 3; j++)
                {

                    connectivity[k] = new Dictionary<int, int>() { { 1, j}, { 2, j + 1} };
                    k += 1;

                }
                connectivity[k] = new Dictionary<int, int>() { { 1, 4 }, { 2, 5 } };
                k += 1;
            //Contact element
            connectivity[k] = new Dictionary<int, int>() { { 1, 3 }, { 2, 4 } };
            return connectivity;
            }

            private static Dictionary<int, bool[]> CreateNodeFAT()
            {
                Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
                for (int i = 1; i <= 5; i++)
                {
                    nodeFAT[i] = new bool[] { true, true, false, false, false, false };
                }
                return nodeFAT;
            }
            private static Dictionary<int, IElementProperties> CreateElementProperties()
            {
                double E = YoungMod;
                double A = area;
                string type = "Bar2D";
                string type2 = "ContactNtN2D";

                Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
                for (int i = 1; i <= 3; i++)
                {
                    elementProperties[i] = new ElementProperties(E, A, type);
                }

                for (int i = 1; i <= 3; i++)
                {
                    elementProperties[i].Density = density;
                }

            elementProperties[4] = new ElementProperties(E, A, type2);
            elementProperties[4].Density = density;
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
            public static Results RunStaticExample()
            {
                #region Structural
                IAssembly elementsAssembly = CreateAssembly();
                elementsAssembly.CreateElementsAssembly();
                elementsAssembly.ActivateBoundaryConditions = true;
                double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();
                int countContactElements = elementsAssembly.CountElementsOfSameType(typeof(ContactNtN2D));
                ShowToGUI.PlotInitialGeometry(elementsAssembly);
                structuralSolution.LinearScheme = new LUFactorization();
                structuralSolution.NonLinearScheme.Tolerance = 1e-5;
                structuralSolution.ActivateNonLinearSolver = true;
                structuralSolution.NonLinearScheme.numberOfLoadSteps = 5;

                double[] externalForces3 = externalForcesStructuralVector;
                foreach (var dof in loadedStructuralDOFs)
                {
                    externalForces3[dof - 1] = externalStructuralLoad;
                }
                double[] reducedExternalForces3 = BoundaryConditionsImposition.ReducedVector(externalForces3, elementsAssembly.BoundedDOFsVector);
                structuralSolution.AssemblyData = elementsAssembly;
                structuralSolution.Solve(reducedExternalForces3);
                double[] solvector3 = structuralSolution.GetSolution();
                elementsAssembly.UpdateDisplacements(solvector3);
                ShowToGUI.PlotFinalGeometry(elementsAssembly);
                double[] fullSolVector3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector3, elementsAssembly.BoundedDOFsVector);
                Dictionary<int, INode> finalNodes = Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullSolVector3);
                double[] xFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item1;
                double[] yFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item2;
                Dictionary<int, double[]> allStepsSolutions = structuralSolution.GetAllStepsSolutions();
                Dictionary<int, double[]> allStepsFullSolutions = new Dictionary<int, double[]>();
                Dictionary<int, Dictionary<int, double[]>> allStepsContactForces = new Dictionary<int, Dictionary<int, double[]>>();
                Dictionary<int, double[]> elementsInternalContactForcesVector;

                for (int i = 1; i <= allStepsSolutions.Count; i++)
                {
                    elementsInternalContactForcesVector = new Dictionary<int, double[]>();
                    elementsAssembly.UpdateDisplacements(allStepsSolutions[i]);
                    elementsInternalContactForcesVector[4] = elementsAssembly.ElementsAssembly[4].CreateInternalGlobalForcesVector();
                    allStepsContactForces[i] = elementsInternalContactForcesVector;
                    string name = "ContactForce" + i.ToString() + ".dat";
                    VectorOperations.PrintVectorToFile(allStepsContactForces.Single(m => m.Key == i).Value.Single(n=>n.Key == 4).Value, @"C:\Users\Public\Documents\" + name);
                }

                for (int i =0; i < allStepsSolutions.Count; i++)
                {
                    allStepsFullSolutions.Add(i+1, BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions.Single(m=>m.Key == i+1).Value, elementsAssembly.BoundedDOFsVector));
                     int j = i + 1;
                    string name = "solution" + j.ToString()+".dat";
                    VectorOperations.PrintVectorToFile(allStepsFullSolutions.Single(m=>m.Key== i + 1).Value, @"C:\Users\Public\Documents\"+ name);
            }
                List<double[]> structuralSolutions = new List<double[]>();
                
                #endregion
                return new Results() { NonlinearSolution = structuralSolutions, SelectedDOF = 2, SolutionType = "Nonlinear" };
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
