using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public static class BendingBeamContact2d
    {

        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        const double gap = 0.01;
        const double thickness = 0.1;
        const int nodesInX1 = 21;
        const int nodesInX2 = 11;
        const int nodesInX3 = 11;
        const int nodesInY1 = 9;
        //const int nodesInY2 = 4;
        const int nodesNumber = 266;
        const int elementsNumber = 70;
        //const int elementsOfType1Number = 40;
        //const int elementsOfType2Number = 30;
        const int contactElements = 8;
        const double xInterv1 = 0.20;
        const double xInterv2 = 0.10;
        const double xInterv3 = 0.20;
        const double yInterv = 0.10;
        const double offset = 0.0;

        //External loads
        const double externalStructuralLoad = -25.0;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;

        const double YoungMod = 1.0 * 1e5;

        const double poissonRatio = 0.25;
        const double density = 8000.0;
        const double area = 1.0;
        const double contactArea = thickness * xInterv2;



        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            boundedDofs.Add(129);
            boundedDofs.Add(130);
            boundedDofs.Add(169);
            boundedDofs.Add(170);
            boundedDofs.Add(299);
            boundedDofs.Add(300);
            boundedDofs.Add(339);
            boundedDofs.Add(340);
            boundedDofs.Add(341);
            boundedDofs.Add(342);
            boundedDofs.Add(363);
            boundedDofs.Add(364);
            boundedDofs.Add(405);
            boundedDofs.Add(406);
            boundedDofs.Add(427);
            boundedDofs.Add(428);
            boundedDofs.Add(469);
            boundedDofs.Add(470);
            boundedDofs.Add(491);
            boundedDofs.Add(492);
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            loadedStructuralDOFs.Add(278);
            externalForcesStructuralVector = new double[nodesNumber * 2];
        }

        private static Dictionary<int, INode> CreateNodes()
        {

            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int k;
            k = 1;
            //First beam
            for(int i = 0; i < 4; i++)
            {
                for (int j = 0; j < nodesInX1; j++)
                {
                    nodes[k] = new Node(j * xInterv2, i * yInterv);
                    k += 1;
                }
                for (int j = 0; j < nodesInX2; j++)
                {
                    nodes[k] = new Node(j * xInterv1, i * yInterv + yInterv / 2.0);
                    k += 1;
                }
            }
            for (int i = 0; i < nodesInX1; i++)
            {
                nodes[k] = new Node(i * xInterv2, yInterv *4);
                k += 1;
            }
            //Second beam
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < nodesInX1; j++)
                {
                    nodes[k] = new Node(j * xInterv2, i * yInterv - gap - (4 - 1) * yInterv);
                    k += 1;
                }
                for (int j = 0; j < nodesInX2; j++)
                {
                    nodes[k] = new Node(j * xInterv1, i * yInterv + yInterv / 2.0 - gap - (4 - 1) * yInterv);
                    k += 1;
                }
            }
            for (int i = 0; i < nodesInX1; i++)
            {
                nodes[k] = new Node(i * xInterv2, yInterv * 3 - gap - (4 - 1) * yInterv);
                k += 1;
            }
            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int k = 1;
            for(int i = 1; i <= 4; i++)
            {
                for(int j = 1; j <= 10; j++)
                {
                    int first = (i - 1) * nodesInX1 + (i - 1) * nodesInX2 + (j - 1) * 2 + 1;
                    int fourth = i * nodesInX1 + (i - 1) * nodesInX2 + j + 1;
                    int fifth = i * nodesInX1 + i* nodesInX2 + (j - 1) * 2 + 3;

                    connectivity[k] = new Dictionary<int, int>() { { 1, first }, { 2, first + 1 }, { 3, first + 2 },
                        { 4, fourth },{ 5, fifth },{ 6, fifth - 1 },{ 7, fifth - 2 },{ 8, fourth - 1 } };
                    k += 1;
                }
            }
            for (int i = 1; i <= 3; i++)
            {
                for (int j = 1; j <= 10; j++)
                {
                    int first = (i - 1) * nodesInX1 + (i - 1) * nodesInX2 + (j - 1) * 2 + 150;
                    int fourth = i * nodesInX1 + (i - 1) * nodesInX2 + j + 150;
                    int fifth = i * nodesInX1 + i * nodesInX2 + (j - 1) * 2 + 152;

                    connectivity[k] = new Dictionary<int, int>() { { 1, first }, { 2, first + 1 }, { 3, first + 2 },
                        { 4, fourth },{ 5, fifth },{ 6, fifth - 1 },{ 7, fifth - 2 },{ 8, fourth - 1 } };
                    k += 1;
                }
            }
            //Contact elements
            int firstMasterNode = 248;
            for (int i = 3; i <= 17; i+=2)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, firstMasterNode }, { 2, firstMasterNode + 1 }, { 3, firstMasterNode + 2 },
                                                             {4, i }, {5, i + 1 }, {6, i + 2 } };
                k += 1;
                firstMasterNode += 2;
            }
            return connectivity;
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
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;

            double A = area;
            double CA = contactArea;

            string type = "Quad8";
            string type3 = "ContactStS2D";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= 40; i++)
            {
                elementProperties[i] = new ElementProperties(E, poissonRatio, A, thickness, density, type);

            }
            for (int i = 41; i <= elementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E, poissonRatio, A, thickness, density, type);
            }
            for(int i =elementsNumber + 1; i<= elementsNumber + contactElements; i++)
            {
                elementProperties[i] = new ElementProperties(E, CA, type3, 10.0, 10, 2, 2);
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
            int countContactElements = elementsAssembly.CountElementsOfSameType(typeof(ContactStS2D));
            ShowToGUI.PlotInitialGeometry(elementsAssembly);
            structuralSolution.LinearScheme = new LUFactorization();
            structuralSolution.NonLinearScheme.Tolerance = 1e-5;
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 30;

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
                for(int j = 1;j<=contactElements; j++)
                {
                    elementsInternalContactForcesVector[elementsNumber + j] = elementsAssembly.ElementsAssembly[elementsNumber + j].CreateInternalGlobalForcesVector();
                }
                allStepsContactForces[i] = elementsInternalContactForcesVector;
                string name = "ContactForces" + i.ToString() + ".dat";
                double[] Vector = new double[contactElements*12];
                int count = 0;
                for (int j = 1; j <= contactElements; j++)
                {
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[0];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[1];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[2];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[3];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[4];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[5];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[6];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[7];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[8];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[9];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[10];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j).Value[11];
                    count += 1;
                }
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
