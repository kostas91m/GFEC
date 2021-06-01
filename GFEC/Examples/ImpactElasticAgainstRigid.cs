using OpenTK.Graphics.ES11;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class ImpactElasticAgainstRigid
    {
        //geometry & initial mesh
        private const double width = 2.0;
        private const double height = 1.0;
        private const double thickness = 0.1;
        private const double xIntervals = 0.1;
        private const double yIntervals = 0.1;
        private const double gap = 0.2;
        private const double angleInDegrees = 45;
        private const double angle = (Math.PI / 180) * angleInDegrees;
        private const double scaleFactor = 0.1;


        //material
        private const double density = 3514;
        private const double density2 = 8000;

        private const double YoungMod = 1.050 * 1e12;
        private const double YoungMod2 = 1 * 1e9;


        private const int nodesInX = (int)(width / xIntervals) + 1;
        private const int nodesInY = (int)(height / yIntervals) + 1;
        private const int nodesNumber = nodesInX * nodesInY + nodesInY * nodesInY;
        private const int elementsNumber = (nodesInX - 1) * (nodesInY - 1) + (nodesInY - 1) * (nodesInY - 1);
        private const int contacts = 2 * nodesInY - 1;

        //external loads & boundary conditions
        private const double externalForce = 0;
        static int[] structuralBoundaryConditions;
        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;
        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int k = 1;
            for (int i = 0; i < nodesInY; i++)
            {
                for (int j = 0; j < nodesInY; j++)
                {
                    nodes[k] = new Node(width/2 + Math.Cos(angle) * j * xIntervals * scaleFactor - Math.Sin(angle) * i * yIntervals * scaleFactor, Math.Cos(angle) * i * yIntervals * scaleFactor + Math.Sin(angle) * j * xIntervals * scaleFactor);// upper body nodes
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
            for (int i = 1 + nodesInY * nodesInY; i <= nodesInX + nodesInY * nodesInY; i++)
            {
                if (i == 1 + nodesInY * nodesInY ||
                    i == nodesInX + nodesInY * nodesInY)
                {
                    boundedDofs.Add(2 * i); //lower side support lower body
                }
                else
                {
                    boundedDofs.Add(2 * i - 1); //lower side support lower body
                    boundedDofs.Add(2 * i); //lower side support lower body
                }
            }
            for (int i = 1 + nodesInY * nodesInY; i <= nodesInX * (nodesInY - 1) + 1 + nodesInY * nodesInY; i += nodesInX)
            {
                boundedDofs.Add(2 * i - 1); //left side support lower body
            }
            for (int i = nodesInX + nodesInY * nodesInY; i <= nodesNumber; i += nodesInX)
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
            for (int j = 0; j < nodesInY - 1; j++)
            {
                for (int i = 1; i < nodesInY; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + j * nodesInY }, { 2, i + 1 + j * nodesInY }, { 3, i + 1 + nodesInY + j * nodesInY }, { 4, i + nodesInY + j * nodesInY } };// upper body
                    k += 1;
                }
            }
            for (int j = 0; j < nodesInY - 1; j++)
            {
                for (int i = 1; i < nodesInX; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + j * nodesInX + nodesInY * nodesInY }, { 2, i + 1 + j * nodesInX + nodesInY * nodesInY }, { 3, i + 1 + nodesInX + j * nodesInX + nodesInY * nodesInY }, { 4, i + nodesInX + j * nodesInX + nodesInY * nodesInY } };//lower body
                    k += 1;
                }
            }
            k = connectivity.Count + 1;
            int upperNode = nodesInY * (nodesInY - 1) + 1;
            for (int i = 1; i < nodesInY; i++)
            {
                int lowerRightNode = nodesNumber - (int)(width / xIntervals) / 2;
                int lowerLeftNode = lowerRightNode - 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };//contacts
                upperNode -= nodesInY;
                k += 1;
            }
            connectivity[k] = new Dictionary<int, int>() { { 1, nodesNumber - (int)(width / xIntervals) / 2 }, { 2, upperNode } };//contacts
            k += 1;
            upperNode = 2;
            for (int i = 2; i <= nodesInY; i++)
            {
                int lowerLeftNode = nodesNumber - (int)(width / xIntervals) / 2;
                int lowerRightNode = lowerLeftNode + 1;
                connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };//contacts
                upperNode +=1;
                k += 1;
            }
            return connectivity;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            double E2 = YoungMod2;
            double A = thickness * yIntervals;
            double contactElemArea1 = thickness * yIntervals * scaleFactor;
            double contactElemArea2 = thickness * xIntervals * scaleFactor;
            string type = "Quad4";
            string type2 = "ContactNtS2D";
            string type3 = "ContactNtN2D";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= (nodesInY - 1) * (nodesInY - 1); i++)
            {
                elementProperties[i] = new ElementProperties(E, A * scaleFactor, type);
            }
            for (int i = 1; i <= (nodesInY - 1) * (nodesInY - 1); i++)
            {
                elementProperties[i].Density = density;
                elementProperties[i].Thickness = thickness;
            }
            for (int i = (nodesInY - 1) * (nodesInY - 1) + 1; i <= elementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(E2, A, type);
            }
            for (int i = (nodesInY - 1) * (nodesInY - 1) + 1; i <= elementsNumber; i++)
            {
                elementProperties[i].Density = density2;
                elementProperties[i].Thickness = thickness;
            }
            for (int i = elementsNumber + 1; i <= elementsNumber + nodesInY; i++)
            {
                if (i == elementsNumber + 1)
                {
                    elementProperties[i] = new ElementProperties(E2, contactElemArea1 / 2, type2);
                    elementProperties[i].Density = density;
                    elementProperties[i].Thickness = thickness;
                }
                else if(i == elementsNumber + nodesInY)
                {
                    elementProperties[i] = new ElementProperties(E2, contactElemArea1 / 2 + contactElemArea2 / 2, type3);
                    elementProperties[i].Density = density;
                    elementProperties[i].Thickness = thickness;
                }
                else
                {
                    elementProperties[i] = new ElementProperties(E2, contactElemArea1, type2);
                    elementProperties[i].Density = density;
                    elementProperties[i].Thickness = thickness;
                }
            }
            for (int i = elementsNumber + nodesInY + 1; i <= elementsNumber + contacts; i++)
            {
                if (i == elementsNumber + contacts)
                {
                    elementProperties[i] = new ElementProperties(E2, contactElemArea2 / 2, type2);
                    elementProperties[i].Density = density;
                    elementProperties[i].Thickness = thickness;
                }
                else
                {
                    elementProperties[i] = new ElementProperties(E2, contactElemArea2, type2);
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
            for (int i = 1; i < 2 * nodesInY * nodesInY - 1; i += 2)
            {
                VelocityVector[i] = -20.0;
            }
            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = BoundaryConditionsImposition.ReducedVector(AccelerationVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialDisplacementVector = BoundaryConditionsImposition.ReducedVector(DisplacementVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialVelocityVector = BoundaryConditionsImposition.ReducedVector(VelocityVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialTime = 0.0;
            //ExplicitSolver newSolver = new ExplicitSolver(0.0534, 100);
            ExplicitSolver newSolver = new ExplicitSolver(0.018, 3600);
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
            double[] fullDynamicSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[3590], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[3591], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[3592], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[3593], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[3594], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol6 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[3595], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol7 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[3596], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol8 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[3597], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol9 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[3598], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol10 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[3599], elementsAssembly.BoundedDOFsVector);
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
