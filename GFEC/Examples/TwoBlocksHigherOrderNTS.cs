using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public static class TwoBlocksHigherOrderNTS
    {

        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        const double length = 1.0;
        const double gap = 0.015;
        const double f = 0.050;
        const double thickness = 0.01;
        const int nodesInX = 7;
        const double nodesInY = 3;

        //External loads
        const double externalStructuralLoad = -10.0;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;

        const double YoungMod = 1.0 * 1e5;
        const double YoungMod2 = 80.0;

        const double poissonRatio = 0.3;
        const double density = 8000.0;
        const double area = 1.0;
        const double contactArea = 0.005;



        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            boundedDofs.Add(1);
            boundedDofs.Add(13);
            boundedDofs.Add(15);
            boundedDofs.Add(27);
            boundedDofs.Add(29);
            boundedDofs.Add(41);

            boundedDofs.Add(43);
            boundedDofs.Add(44);
            boundedDofs.Add(46);
            boundedDofs.Add(48);

            boundedDofs.Add(59);
            boundedDofs.Add(60);
            boundedDofs.Add(61);
            boundedDofs.Add(62);
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            loadedStructuralDOFs.Add(34);
            loadedStructuralDOFs.Add(36);
            loadedStructuralDOFs.Add(38);
            externalForcesStructuralVector = new double[31 * 2];
        }

        private static Dictionary<int, INode> CreateNodes()
        {

            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int k;
            k = 1;
            //First block
            nodes[k] = new Node(0.0, 0.0);
            k += 1;
            nodes[k] = new Node(length/2.0, 0.0);
            k += 1;
            nodes[k] = new Node(length, 0.0);
            k += 1;
            nodes[k] = new Node(3.0 * length / 2.0, 0.0);
            k += 1;
            nodes[k] = new Node(2.0 * length, 0.0);
            k += 1;
            nodes[k] = new Node(5.0 * length / 2.0, 0.0);
            k += 1;
            nodes[k] = new Node(3.0 * length, 0.0);
            k += 1;
            nodes[k] = new Node(0.0, length / 2.0);
            k += 1;
            nodes[k] = new Node(length / 2.0, length / 2.0);
            k += 1;
            nodes[k] = new Node(length, length / 2.0);
            k += 1;
            nodes[k] = new Node(3.0 * length / 2.0, length / 2.0);
            k += 1;
            nodes[k] = new Node(2.0 * length, length / 2.0);
            k += 1;
            nodes[k] = new Node(5.0 * length / 2.0, length / 2.0);
            k += 1;
            nodes[k] = new Node(3.0 * length, length / 2.0);
            k += 1;
            nodes[k] = new Node(0.0, length);
            k += 1;
            nodes[k] = new Node(length / 2.0, length);
            k += 1;
            nodes[k] = new Node(length, length);
            k += 1;
            nodes[k] = new Node(3.0 * length / 2.0, length);
            k += 1;
            nodes[k] = new Node(2.0 * length, length);
            k += 1;
            nodes[k] = new Node(5.0 * length / 2.0, length);
            k += 1;
            nodes[k] = new Node(3.0 * length, length);
            k += 1;
            //Second block
            nodes[k] = new Node(0.50 * length, -length - gap - f);
            k += 1;
            nodes[k] = new Node(1.5 * length, -length - gap - f);
            k += 1;
            nodes[k] = new Node(2.5 * length, -length - gap - f);
            k += 1;
            nodes[k] = new Node(2.5 * length, -length/2.0 - gap - f);
            k += 1;
            nodes[k] = new Node(2.5 * length, - gap - f);
            k += 1;
            nodes[k] = new Node(1.5 * length, -gap);
            k += 1;
            nodes[k] = new Node(0.50 * length, -gap - f);
            k += 1;
            nodes[k] = new Node(0.50 * length, -length / 2.0 - gap - f);
            k += 1;
            //Bars
            nodes[k] = new Node(0.0, 2 * length);
            k += 1;
            nodes[k] = new Node(3 * length, 2 * length);
            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int k = 1;
            for(int j = 1; j <= nodesInY - 1; j++)
            {
                for(int i = 1; i<=nodesInX - 1; i++)
                {
                    int count = (j - 1) * nodesInX;
                    connectivity[k] = new Dictionary<int, int>() { { 1, count + i }, { 2, count + i + 1 }, { 3, count + i + 1 + nodesInX }, { 4, count + i + nodesInX } };
                    k += 1;
                }
            }
            connectivity[k] = new Dictionary<int, int>() { { 1, 22 }, { 2, 23 }, { 3, 24 }, { 4, 25 },
                                                         { 5, 26 }, { 6, 27 }, { 7, 28 }, { 8, 29 }};
            k += 1;
            //bar elements
            connectivity[k] = new Dictionary<int, int>() { { 1, 15 }, { 2, 30 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 21 }, { 2, 31 } };
            k += 1;
            //Contact elements
            for(int i = 2; i<=6; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, 28 }, { 2, 27 },
                                                         { 3, 26 }, { 4, i }};
                k += 1;
            }
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= 31; i++)
            {
                nodeFAT[i] = new bool[] { true, true, false, false, false, false };
            }
            return nodeFAT;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            double E2 = YoungMod2;

            double A = area;
            double CA = contactArea;

            string type = "Quad4";
            string type2 = "Quad8";
            string type3 = "SecondOrderContactNtS2D";
            string type4 = "Bar2D";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= 12; i++)
            {
                elementProperties[i] = new ElementProperties(E, poissonRatio, A, thickness, density, type);

            }
            elementProperties[13] = new ElementProperties(E, poissonRatio, A, thickness, density, type2);
            for (int i = 14; i <= 15; i++)
            {
                elementProperties[i] = new ElementProperties(E2, A, type4);
                elementProperties[i].Density = density;
            }
            for(int i = 16; i<=20; i++)
            {
                elementProperties[i] = new ElementProperties(E, CA, type3);
                elementProperties[i].Density = density;
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
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 50;

            double[] externalForces3 = externalForcesStructuralVector;
            foreach (var dof in loadedStructuralDOFs)
            {
                externalForces3[dof - 1] = externalStructuralLoad;
            }
            double[] reducedExternalForces3 = BoundaryConditionsImposition.ReducedVector(externalForces3, elementsAssembly.BoundedDOFsVector);
            structuralSolution.AssemblyData = elementsAssembly;
            structuralSolution.Solve(reducedExternalForces3);
            double[] solvector3 = structuralSolution.GetSolution();
            Dictionary<int, double[]> allStepsSolutions = structuralSolution.GetAllStepsSolutions();
            Dictionary<int, List<double[]>> gPointsStress = new Dictionary<int, List<double[]>>();
            Dictionary<int, List<double[]>> gPointsStrain = new Dictionary<int, List<double[]>>();
            Dictionary<int, List<double[]>> gPoints = new Dictionary<int, List<double[]>>();
            Dictionary<int, List<double[]>> nodalStress = new Dictionary<int, List<double[]>>();
            Dictionary<int, List<double[]>> nodalStrain = new Dictionary<int, List<double[]>>();
            for (int i = 1; i <= allStepsSolutions.Count; i++)
            {
                string name = "NodalCoordinates" + i.ToString() + ".dat";
                ExportToFile.ExportUpdatedNodalCoordinates(elementsAssembly, BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions.Single(m => m.Key == i).Value, elementsAssembly.BoundedDOFsVector), name);
                gPointsStress = elementsAssembly.GetElementsStresses(allStepsSolutions[i]);
                gPointsStrain = elementsAssembly.GetElementsStains(allStepsSolutions[i]);
                gPoints = elementsAssembly.GetElementsGaussPoints(allStepsSolutions[i]);
                nodalStress = elementsAssembly.GetElementsNodesStresses(allStepsSolutions[i]);
                nodalStrain = elementsAssembly.GetElementsNodesStains(allStepsSolutions[i]);
                string name1 = "GPointsStress" + i.ToString() + ".dat";
                string name2 = "GPointsStrain" + i.ToString() + ".dat";
                string name3 = "GPointsCoordinates" + i.ToString() + ".dat";
                string name4 = "NodalStress" + i.ToString() + ".dat";
                string name5 = "NodalStrain" + i.ToString() + ".dat";

                VectorOperations.PrintDictionaryofListsofVectorsToFile(gPointsStress, @"C:\Users\Public\Documents\" + name1);
                VectorOperations.PrintDictionaryofListsofVectorsToFile(gPointsStrain, @"C:\Users\Public\Documents\" + name2);
                VectorOperations.PrintDictionaryofListsofVectorsToFile(gPoints, @"C:\Users\Public\Documents\" + name3);
                VectorOperations.PrintDictionaryofListsofVectorsToFile(nodalStress, @"C:\Users\Public\Documents\" + name4);
                VectorOperations.PrintDictionaryofListsofVectorsToFile(nodalStrain, @"C:\Users\Public\Documents\" + name5);
            }
            elementsAssembly.UpdateDisplacements(solvector3);
            ShowToGUI.PlotFinalGeometry(elementsAssembly);
            double[] fullSolVector3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector3, elementsAssembly.BoundedDOFsVector);
            Dictionary<int, INode> finalNodes = Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullSolVector3);
            double[] xFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item1;
            double[] yFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item2;
            Dictionary<int, double[]> allStepsFullSolutions = new Dictionary<int, double[]>();
            Dictionary<int, Dictionary<int, double[]>> allStepsContactForces = new Dictionary<int, Dictionary<int, double[]>>();
            Dictionary<int, double[]> elementsInternalContactForcesVector;
            for (int i = 1; i <= allStepsSolutions.Count; i++)
            {
                elementsInternalContactForcesVector = new Dictionary<int, double[]>();
                elementsAssembly.UpdateDisplacements(allStepsSolutions[i]);
                elementsInternalContactForcesVector[16] = elementsAssembly.ElementsAssembly[16].CreateInternalGlobalForcesVector();
                elementsInternalContactForcesVector[17] = elementsAssembly.ElementsAssembly[17].CreateInternalGlobalForcesVector();
                elementsInternalContactForcesVector[18] = elementsAssembly.ElementsAssembly[18].CreateInternalGlobalForcesVector();
                elementsInternalContactForcesVector[19] = elementsAssembly.ElementsAssembly[19].CreateInternalGlobalForcesVector();
                elementsInternalContactForcesVector[20] = elementsAssembly.ElementsAssembly[20].CreateInternalGlobalForcesVector();
                allStepsContactForces[i] = elementsInternalContactForcesVector;
                string name = "ContactForces" + i.ToString() + ".dat";
                double[] Vector = new double[40];
                int count = 0;
                int count2 = 16;
                for(int j =0; j< 5; j++)
                {
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == count2).Value[0];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == count2).Value[1];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == count2).Value[2];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == count2).Value[3];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == count2).Value[4];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == count2).Value[5];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == count2).Value[6];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == count2).Value[7];
                    count += 1;
                    count2 += 1;
                }
                //Vector[8] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 8).Value[0];
                //Vector[9] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 8).Value[1];
                //Vector[10] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 8).Value[2];
                //Vector[11] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 8).Value[3];
                //Vector[12] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 8).Value[4];
                //Vector[13] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 8).Value[5];
                //Vector[14] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 8).Value[6];
                //Vector[15] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == 8).Value[7];
                VectorOperations.PrintVectorToFile(Vector, @"C:\Users\Public\Documents\" + name);
            }

            for (int i = 0; i < allStepsSolutions.Count; i++)
            {
                allStepsFullSolutions.Add(i + 1, BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions.Single(m => m.Key == i + 1).Value, elementsAssembly.BoundedDOFsVector));
                int j = i + 1;
                string name = "solution" + j.ToString() + ".dat";
                VectorOperations.PrintVectorToFile(allStepsFullSolutions.Single(m => m.Key == i + 1).Value, @"C:\Users\Public\Documents\" + name);
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
