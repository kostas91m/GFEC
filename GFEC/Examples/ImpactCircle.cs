using OpenTK.Graphics.ES11;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class ImpactCircle
    {
        //geometry & initial mesh
        private const double width = 2.0;
        private const double height = 1.0;
        private const double thickness = 0.1;
        private const double xIntervals = 0.1;
        private const double yIntervals = 0.1;
        //private const double angleInDegrees = 45;
        //private const double angle = (Math.PI / 180) * angleInDegrees;
        private const double Radius = 1.0;
        private const double scaleFactor = 0.33;
        private const double gap = Radius * scaleFactor + 0.02;



        //material
        private const double density = 3514;
        private const double density2 = 8000;

        private const double YoungMod = 1.050 * 1e12;
        private const double YoungMod2 = 2 * 1e9;
        private const double PoissonRatio = 0.33;


        private const int nodesInX = (int)(width / xIntervals) + 1;
        private const int nodesInY = (int)(height / yIntervals) + 1;
        private const int nodesInPerimeter = 36;

        private const int nodesNumber = nodesInX * nodesInY + nodesInPerimeter * 4 + 1 + 12;
        private const int elementsNumber = (nodesInX - 1) * (nodesInY - 1) + 3 * nodesInPerimeter + 60;
        private const int contacts = 11;

        //external loads & boundary conditions
        private const double externalForce = 0;
        static int[] structuralBoundaryConditions;
        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;
        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            nodes[1] = new Node(0, 0);// center
            int k = 2;
            double theta = 0;
            for (int i = 0; i < 12; i++)
            {
                double cos = Math.Cos(theta);
                double sin = Math.Sin(theta);
                nodes[k] = new Node((0.50) * Radius * scaleFactor * cos, (0.50) * Radius * scaleFactor * sin);// upper body nodes
                theta += 30 * (Math.PI / 180);
                k += 1;
            }
            double theta5 = 0;
            for (int i = 0; i < 12; i++)
            {
                double cos5 = Math.Cos(theta5);
                double sin5 = Math.Sin(theta5);
                nodes[k] = new Node((0.666667) * Radius * scaleFactor * cos5, (0.666667) * Radius * scaleFactor * sin5);// upper body nodes
                theta5 +=  10 * (Math.PI / 180);
                cos5 = Math.Cos(theta5);
                sin5 = Math.Sin(theta5);
                k += 1;
                nodes[k] = new Node((0.6467335) * Radius * scaleFactor * cos5, (0.6467335) * Radius * scaleFactor * sin5);// upper body nodes
                theta5 +=  10 * (Math.PI / 180);
                cos5 = Math.Cos(theta5);
                sin5 = Math.Sin(theta5);
                k += 1;
                nodes[k] = new Node((0.6467335) * Radius * scaleFactor * cos5, (0.6467335) * Radius * scaleFactor * sin5);// upper body nodes
                k += 1;
                theta5 += 10 * (Math.PI / 180);
            }
            double theta2 = 0;
            for (int i = 0; i < nodesInPerimeter; i++)
            {
                double cos2 = Math.Cos(theta2);
                double sin2 = Math.Sin(theta2);
                nodes[k] = new Node(0.75 * Radius * scaleFactor * cos2, 0.75 * Radius * scaleFactor * sin2);// upper body nodes
                theta2 +=  10 * (Math.PI / 180);
                k += 1;
            }
            double theta3 = 0;
            for (int i = 0; i < nodesInPerimeter; i++)
            {
                double cos3 = Math.Cos(theta3);
                double sin3 = Math.Sin(theta3);
                nodes[k] = new Node((0.88) * Radius * scaleFactor * cos3, (0.88) * Radius * scaleFactor * sin3);// upper body nodes
                k += 1;
                theta3 +=  10 * (Math.PI / 180);
            }
            double theta4 = 0;
            for (int i = 0; i < nodesInPerimeter; i++)
            {
                double cos4 = Math.Cos(theta4);
                double sin4 = Math.Sin(theta4);
                nodes[k] = new Node(Radius * scaleFactor * cos4, Radius * scaleFactor * sin4);// upper body nodes
                k += 1;
                theta4 += 10 * (Math.PI / 180);
            }
            for (int i = 0; i < nodesInY; i++)
            {
                for (int j = 0; j < nodesInX; j++)
                {
                    nodes[k] = new Node(j * xIntervals - width/2, i * yIntervals - gap - height);// lower body nodes
                    k += 1;
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
            //for (int i = 1; i <= nodesInX * (nodesInY - 1) + 1; i += nodesInX)
            //{
            //    boundedDofs.Add(2 * i - 1); //left side support upper body
            //}
            //for (int i = nodesInX; i <= nodesInX * nodesInY; i += nodesInX)
            //{
            //    boundedDofs.Add(2 * i - 1); //right side support upper body
            //}
            for (int i = 1 + nodesInPerimeter * 4 + 13; i <= nodesInX + nodesInPerimeter * 4 + 13; i++)
            {
                if (i == 1 + nodesInPerimeter * 4 + 13 ||
                    i == nodesInX + nodesInPerimeter * 4 + 13)
                {
                    boundedDofs.Add(2 * i); //lower side support lower body
                }
                else
                {
                    boundedDofs.Add(2 * i - 1); //lower side support lower body
                    boundedDofs.Add(2 * i); //lower side support lower body
                }
            }
            for (int i = 1 + nodesInPerimeter * 4 + 13; i <= nodesInX * (nodesInY - 1) + 1 + nodesInPerimeter * 4 + 13; i += nodesInX)
            {
                boundedDofs.Add(2 * i - 1); //left side support lower body
            }
            for (int i = nodesInX + nodesInPerimeter * 4 + 13; i <= nodesNumber; i += nodesInX)
            {
                boundedDofs.Add(2 * i - 1); //right side support lower body
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
            //externalForcesStructuralVector = new double[nodesInY * nodesInY];

            return externalForcesStructuralVector;
        }
        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int k = 1;
            int second = 2;
            for (int j = 0; j < 12; j++)
            {
                int first = 1;
                int third = second + 1;
                if (j < 11)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, first }, { 2, second }, { 3, third } };// upper body
                    second = third;
                }
                else
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, first }, { 2, second }, { 3, 2 } };// upper body
                }

                k += 1;
            }
            for (int j = 0; j < 12; j++)
            {
                int first2 = j + 2;
                if (j < 11)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, first2 }, { 2, first2 + 12 + j * 2 }, { 3, first2 + 13 + j * 2 } };// upper body
                    k += 1;
                    connectivity[k] = new Dictionary<int, int>() { { 1, first2 }, { 2, first2 + 13 + j * 2 }, { 3, first2 + 1 } };// upper body
                    k += 1;
                    connectivity[k] = new Dictionary<int, int>() { { 1, first2 + 1 }, { 2, first2 + 13 + j * 2 }, { 3, first2 + 14 + j * 2 } };// upper body
                    k += 1;
                    connectivity[k] = new Dictionary<int, int>() { { 1, first2 + 1 }, { 2, first2 + 14 + j * 2 }, { 3, first2 + 15 + j * 2 } };// upper body
                    k += 1;
                }
                else
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, first2 }, { 2, first2 + 12 + j * 2 }, { 3, first2 + 13 + j * 2 } };// upper body
                    k += 1;
                    connectivity[k] = new Dictionary<int, int>() { { 1, first2 }, { 2, first2 + 13 + j * 2 }, { 3, 2 } };// upper body
                    k += 1;
                    connectivity[k] = new Dictionary<int, int>() { { 1, 2 }, { 2, first2 + 13 + j * 2 }, { 3, first2 + 14 + j * 2 } };// upper body
                    k += 1;
                    connectivity[k] = new Dictionary<int, int>() { { 1, 2 }, { 2, first2 + 14 + j * 2 }, { 3, 14 } };// upper body
                    k += 1;
                }
            }
            for (int j = 1; j < nodesInPerimeter; j++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, j + 13 }, { 2, j + 13 + nodesInPerimeter }, { 3, j + 2 + nodesInPerimeter + 12 }, { 4, j + 2 + 12 } };// upper body
                k += 1;
            }
            connectivity[k] = new Dictionary<int, int>() { { 1, 49 }, { 2, 49 + nodesInPerimeter }, { 3, 2 + nodesInPerimeter + 12 }, { 4, 2 + 12 } };// upper body
            k += 1;
            for (int j = 1; j < nodesInPerimeter; j++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, j + 1 + nodesInPerimeter + 12 }, { 2, j + 1 + 2 * nodesInPerimeter + 12 }, { 3, j + 2 + 2 * nodesInPerimeter + 12 }, { 4, j + 2 + nodesInPerimeter + 12 } };// upper body
                k += 1;
            }
            connectivity[k] = new Dictionary<int, int>() { { 1, 37 + nodesInPerimeter + 12 }, { 2, 37 + 2 * nodesInPerimeter + 12 }, { 3, 2 + 2 * nodesInPerimeter + 12 }, { 4, 2 + nodesInPerimeter + 12 } };// upper body
            k += 1;
            for (int j = 1; j < nodesInPerimeter; j++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, j + 1 + 2 * nodesInPerimeter + 12 }, { 2, j + 1 + 3 * nodesInPerimeter + 12 }, { 3, j + 2 + 3 * nodesInPerimeter + 12 }, { 4, j + 2 + 2 * nodesInPerimeter + 12 } };// upper body
                k += 1;
            }
            connectivity[k] = new Dictionary<int, int>() { { 1, 37 + 2 * nodesInPerimeter + 12 }, { 2, 37 + 3 * nodesInPerimeter + 12 }, { 3, 2 + 3 * nodesInPerimeter + 12 }, { 4, 2 + 2 * nodesInPerimeter + 12 } };// upper body
            k += 1;
            for (int j = 0; j < nodesInY - 1; j++)
            {
                for (int i = 1; i < nodesInX; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + j * nodesInX + nodesInPerimeter * 4 + 1 + 12 }, { 2, i + 1 + j * nodesInX + nodesInPerimeter * 4 + 1 + 12 }, { 3, i + 1 + nodesInX + j * nodesInX + nodesInPerimeter * 4 + 1 + 12 }, { 4, i + nodesInX + j * nodesInX + nodesInPerimeter * 4 + 1 + 12 } };//lower body
                    k += 1;
                }
            }
            k = connectivity.Count + 1;
            int centerLowerNode = nodesNumber - (int)(width / xIntervals) / 2;
            int centerTopNode = 149;
            connectivity[k] = new Dictionary<int, int>() { { 1, centerLowerNode }, { 2, centerTopNode } };//contacts
            k += 1;
            int upperNode = centerTopNode - 1;
            int lowerLeftNode = centerLowerNode - 1;
            int lowerRightNode = lowerLeftNode + 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };//contacts
            k += 1;
            upperNode = centerTopNode + 1;
            lowerLeftNode = centerLowerNode;
            lowerRightNode = lowerLeftNode + 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };//contacts
            k += 1;
            upperNode = centerTopNode - 2;
            lowerLeftNode = centerLowerNode - 2;
            lowerRightNode = lowerLeftNode + 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };//contacts
            k += 1;
            upperNode = centerTopNode + 2;
            lowerLeftNode = centerLowerNode + 1;
            lowerRightNode = lowerLeftNode + 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };//contacts
            k += 1;
            upperNode = centerTopNode - 3;
            lowerLeftNode = centerLowerNode - 2;
            lowerRightNode = lowerLeftNode + 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };//contacts
            k += 1;
            upperNode = centerTopNode + 3;
            lowerLeftNode = centerLowerNode + 1;
            lowerRightNode = lowerLeftNode + 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };//contacts
            k += 1;
            upperNode = centerTopNode - 4;
            lowerLeftNode = centerLowerNode - 3;
            lowerRightNode = lowerLeftNode + 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };//contacts
            k += 1;
            upperNode = centerTopNode + 4;
            lowerLeftNode = centerLowerNode + 2;
            lowerRightNode = lowerLeftNode + 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };//contacts
            k += 1;
            upperNode = centerTopNode - 5;
            lowerLeftNode = centerLowerNode - 3;
            lowerRightNode = lowerLeftNode + 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };//contacts
            k += 1;
            upperNode = centerTopNode + 5;
            lowerLeftNode = centerLowerNode + 2;
            lowerRightNode = lowerLeftNode + 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };//contacts
            return connectivity;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            double E2 = YoungMod2;
            double A = thickness * yIntervals;
            double contactElemArea = thickness * 0.17431 * Radius * scaleFactor;
            string type = "Quad4";
            string type2 = "Triangle3";
            //string type2 = "ContactNtS2D";
            string type3 = "ContactNtN2D";
            string type4 = "ContactNtS2D";


            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= 12; i++)
            {
                elementProperties[i] = new ElementProperties(E, PoissonRatio, A, thickness, density, type2);
            }
            for (int i = 13; i <= 60; i++)
            {
                elementProperties[i] = new ElementProperties(E, PoissonRatio, A, thickness, density, type2);
            }
            for (int i = 61; i <= 3 * nodesInPerimeter + 60; i++)
            {
                elementProperties[i] = new ElementProperties(E, PoissonRatio, A, thickness, density2, type);
            }
            for (int i = 3 * nodesInPerimeter + 60 + 1; i <= elementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E2, PoissonRatio, A, thickness, density, type);
            }
            elementProperties[elementsNumber + 1] = new ElementProperties(E2, contactElemArea, type3);
            elementProperties[elementsNumber + 1].Density = density;
            elementProperties[elementsNumber + 1].Thickness = thickness;
            for (int i = elementsNumber + 2; i <= elementsNumber + contacts; i++)
            {
                elementProperties[i] = new ElementProperties(E2, contactElemArea, type4);
                elementProperties[i].Density = density;
                elementProperties[i].Thickness = thickness;
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
            for (int i = 1; i <= 2 * (nodesInPerimeter * 4 + 13) - 1; i += 2)
            {
                VelocityVector[i] = -30.0;
            }
            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = BoundaryConditionsImposition.ReducedVector(AccelerationVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialDisplacementVector = BoundaryConditionsImposition.ReducedVector(DisplacementVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialVelocityVector = BoundaryConditionsImposition.ReducedVector(VelocityVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialTime = 0.0;
            //ExplicitSolver newSolver = new ExplicitSolver(0.0534, 100);
            ExplicitSolver newSolver = new ExplicitSolver(0.002906, 6000);
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
            double[] fullDynamicSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[510], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[520], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[530], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[540], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[550], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol6 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[560], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol7 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[570], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol8 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[580], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol9 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[590], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol10 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[599], elementsAssembly.BoundedDOFsVector);
            VectorOperations.PrintVectorToFile(fullDynamicSol1, @"C:\Users\Public\Documents\Results1.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol2, @"C:\Users\Public\Documents\Results2.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol3, @"C:\Users\Public\Documents\Results3.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol4, @"C:\Users\Public\Documents\Results4.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol5, @"C:\Users\Public\Documents\Results5.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol6, @"C:\Users\Public\Documents\Results6.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol7, @"C:\Users\Public\Documents\Results7.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol8, @"C:\Users\Public\Documents\Results8.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol9, @"C:\Users\Public\Documents\Results9.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol10, @"C:\Users\Public\Documents\Results10.dat");
            //newSolver.PrintExplicitSolution();
            Results finalResults = new Results() { DynamicSolution = newSolver.explicitSolution, TimeSteps = newSolver.TimeAtEachStep, SelectedDOF = 1, SelectedInterval = 1, SolutionType = "Dynamic" };
            return finalResults;
        }
    }
}
