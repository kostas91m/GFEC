using OpenTK.Graphics.ES11;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class ImpactCircle2
    {
        //geometry & initial mesh
        private const double width = 4;
        private const double height = 2;
        private const double thickness = 0.1;
        private const double xIntervals = 0.2;
        private const double yIntervals = 0.2;
        private const double gap = 1.01;
        private const int circles = 3;
        private const int steps = 36;
        private const double Radius = 1.00;
        private const double radius = 0.70;
        private const double RD = (Radius - radius) / (circles - 1);
        private const double angle = (Math.PI / 180) * (360 / steps);


        //material
        private const double density = 8000;
        private const double massScaleFactor = 1.0;
        private const double density2 = 8000;

        private const double YoungMod = 300.0 * 1e9;
        private const double YoungMod2 = 5.0 * 1e9;


        private const int nodesInX = (int)(width / xIntervals) + 1;
        private const int nodesInY = (int)(height / yIntervals) + 1;
        //private const int nodesNumber = nodesInX * nodesInY + circles * steps;
        private const int nodesNumber = nodesInX * nodesInY + (3*circles - 1) * steps;
        private const int elementsNumber = (nodesInX - 1) * (nodesInY - 1) + (circles - 1) * (steps);
        private const int contacts = 13;

        //external loads & boundary conditions
        private const double externalForce = 0;
        static int[] structuralBoundaryConditions;
        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;
        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int k = 1;
            for (int i = circles; i > 0; i--)
            {
                for (int j = 0; j < steps; j++)
                {
                    nodes[k] = new Node((radius + RD * (i - 1)) * Math.Cos(j * (angle)), ((radius + RD * (i - 1)) * Math.Sin(j * (angle))));// upper body nodes
                    k += 1;
                }
            }
            double initialAngle = angle / 2.0;
            for (int i = circles; i > 0; i--)
            {
                for (int j = 0; j < steps; j++)
                {
                    nodes[k] = new Node((radius + RD * (i - 1)) * Math.Cos(j * (angle) + initialAngle), ((radius + RD * (i - 1)) * Math.Sin(j * (angle) + initialAngle)));// upper body nodes
                    k += 1;
                }
            }
            double radialDisplacement = RD / 2.0;
            for (int i = circles; i > 1; i--)
            {
                for (int j = 0; j < steps; j++)
                {
                    nodes[k] = new Node((radius + RD * (i - 1) - radialDisplacement) * Math.Cos(j * (angle)), ((radius + RD * (i - 1) - radialDisplacement) * Math.Sin(j * (angle))));// upper body nodes
                    k += 1;
                }
            }
            for (int i = 0; i < nodesInY; i++)
            {
                for (int j = 0; j < nodesInX; j++)
                {
                    nodes[k] = new Node(j * xIntervals - (width / 2), i * yIntervals - gap - height);// lower body nodes
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
            for (int i = 1 + (3 * circles - 1) * steps; i <= nodesInX + (3 * circles - 1) * steps; i++)
            {
                if (i == 1 + (3 * circles - 1) * steps ||
                    i == nodesInX + (3 * circles - 1) * steps)
                {
                    boundedDofs.Add(2 * i); //lower side support lower body
                }
                else
                {
                    boundedDofs.Add(2 * i - 1); //lower side support lower body
                    boundedDofs.Add(2 * i); //lower side support lower body
                }
            }
            for (int i = 1 + (3 * circles - 1) * steps; i <= nodesInX * (nodesInY - 1) + 1 + (3 * circles - 1) * steps; i += nodesInX)
            {
                boundedDofs.Add(2 * i - 1); //left side support lower body
            }
            for (int i = nodesInX + (3 * circles - 1) * steps; i <= nodesNumber; i += nodesInX)
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
            for (int j = 0; j < circles - 1; j++)
            {
                int ml = 1 + j * steps;
                int mu = ml + steps;
                int mr1 = 2 + j* steps;
                int mr2 = 1 + j*steps;
                for (int i = 1; i < steps; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + j * steps }, { 2, ml + circles *steps }, { 3, i + 1 + j * steps }, { 4, mr1 + 2*circles *steps },
                        { 5, i + 1 + steps + j * steps }, { 6, mu + circles *steps }, { 7, i + steps + j * steps }, { 8,2*circles *steps + mr2 } };// upper body
                    k += 1;
                    ml += 1;
                    mu += 1;
                    mr1 += 1;
                    mr2 += 1;
                }
            }
            for (int i = 1; i < circles; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i * steps }, { 2, steps + circles *steps +(i-1)*steps }, { 3, (i - 1) * steps + 1}, { 4, 1 +(i-1)*steps + 2*circles *steps },
                    { 5, i * steps + 1  }, { 6, 2*steps + circles *steps +(i-1)*steps }, { 7, (i + 1) * steps }, { 8, steps +(i-1)*steps + 2*circles *steps } };// close circle upper body

                k += 1;
            }
            for (int j = 0; j < nodesInY - 1; j++)
            {
                for (int i = 1; i < nodesInX; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + j * nodesInX + (3 * circles - 1) * steps }, { 2, i + 1 + j * nodesInX + (3 * circles - 1) * steps }, { 3, i + 1 + nodesInX + j * nodesInX + (3 * circles - 1) * steps }, { 4, i + nodesInX + j * nodesInX + (3 * circles - 1) * steps } };//lower body
                    k += 1;
                }
            }
            k = connectivity.Count + 1;

            int lowerRightNode = nodesNumber + 1 - (int)(width / xIntervals) / 2;
            int lowerLeftNode = lowerRightNode - 2;
            int lowerCentralNode = lowerRightNode - 1;

            connectivity[k] = new Dictionary<int, int>() { { 1, lowerCentralNode }, { 2, lowerRightNode }, { 3, 3*steps / 4 + 2 } };//contacts
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerRightNode }, { 2, lowerRightNode + 1 }, { 3, 3*steps / 4 + 3 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerRightNode + 1}, { 2, lowerRightNode + 2 }, { 3, 3*steps / 4 + 4 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerCentralNode }, { 3, 3*steps / 4} };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode - 1 }, { 2, lowerCentralNode - 1 }, { 3, 3*steps / 4 - 1 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode - 2 }, { 2, lowerCentralNode - 2 }, { 3, 3*steps / 4 - 2 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerCentralNode }, { 2, 3*steps / 4 + 1 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode - 2 }, { 2, lowerLeftNode - 1 }, { 3, circles * steps + 25 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode - 1 }, { 2, lowerLeftNode }, { 3, circles * steps + 26 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerCentralNode }, { 3, circles * steps + 27 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerCentralNode}, { 2, lowerRightNode}, { 3, circles * steps + 28 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerRightNode }, { 2, lowerRightNode + 1 }, { 3, circles * steps + 29 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, lowerRightNode + 1 }, { 2, lowerRightNode + 2 }, { 3, circles * steps + 30 } };
            return connectivity;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            double E2 = YoungMod2;
            double A = thickness * yIntervals;
            double CircularA = thickness * 2 * Math.PI * Radius / steps;
            //double A1 = thickness * (Math.PI * Radius * Radius - Math.PI * (Radius - ((Radius - radius) / (circles - 1))) * (Radius - ((Radius - radius) / (circles - 1))))/ (steps - 1);
            //double A2 = thickness * (Math.PI * (Radius- RD) * (Radius - RD) - Math.PI * ((Radius - RD) - ((Radius - radius) / (circles - 1))) * ((Radius - RD) - ((Radius - radius) / (circles - 1)))) / (steps - 1);
            //double A3 = thickness * (Math.PI * (Radius - 2 * RD) * (Radius - 2 *RD) - Math.PI * ((Radius - 2 * RD) - ((Radius - radius) / (circles - 1))) * ((Radius - 2 * RD) - ((Radius - radius) / (circles - 1)))) / (steps - 1);
            string type = "Quad8";
            string type2 = "ContactNtS2D";
            string type3 = "ContactNtN2D";
            string type4 = "Quad4";


            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();

            //for (int j = 1; j < circles; j++)                                      // upper body element properties
            //{
            //double CircularSectorArea = thickness * (Math.PI * (Radius - (j - 1) * RD) * (Radius - (j - 1) * RD) - Math.PI * (Radius - j * RD) * (Radius - j * RD)) / steps;

            //int m = (j - 1) * steps + 1;
            for (int i = 1; i <= (circles - 1) * steps; i++)
            {
                elementProperties[i] = new ElementProperties(E, CircularA, type);
                elementProperties[i].Density = density * massScaleFactor;
                elementProperties[i].Thickness = thickness;
            }
            //}                                                                   // end of upper body properties


            for (int i = (circles - 1) * steps + 1; i <= elementsNumber; i++) //Lower Body elements properties
            {
                elementProperties[i] = new ElementProperties(E2, A, type4);
                elementProperties[i].Density = density2;
                elementProperties[i].Thickness = thickness;

            }
            for (int i = elementsNumber + 1; i <= elementsNumber + contacts; i++) //Contacts elements properties
            {
                if (i == elementsNumber + 7)
                {
                    elementProperties[i] = new ElementProperties(E2, CircularA, type3);
                    elementProperties[i].Density = density;
                    elementProperties[i].Thickness = thickness;
                }
                else
                {
                    elementProperties[i] = new ElementProperties(E2, CircularA, type2);
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
            for (int i = 1; i <= 2 * (3 * circles - 1) * steps - 1; i += 2)
            {
                VelocityVector[i] = -50.0;
            }
            InitialConditions initialValues = new InitialConditions();
            initialValues.InitialAccelerationVector = BoundaryConditionsImposition.ReducedVector(AccelerationVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialDisplacementVector = BoundaryConditionsImposition.ReducedVector(DisplacementVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialVelocityVector = BoundaryConditionsImposition.ReducedVector(VelocityVector, elementsAssembly.BoundedDOFsVector);
            initialValues.InitialTime = 0.0;
            //ExplicitSolver newSolver = new ExplicitSolver(0.0035, 3500);
            ExplicitSolver newSolver = new ExplicitSolver(0.003, 300);
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
            newSolver.SolveNewmark();
            Tuple<Dictionary<int, double[]>, Dictionary<int, double>> solvectors = newSolver.GetResults();
            //int max = solvectors.Item1.OrderByDescending(m => m.Key).FirstOrDefault().Key;
            //elementsAssembly.UpdateDisplacements(solvectors.Item1.Single(m => m.Key == max).Value);
            Dictionary<int, double[]> allStepsSolutions = solvectors.Item1;
            Dictionary<int, List<double[]>> stress = new Dictionary<int, List<double[]>>();
            Dictionary<int, List<double[]>> strain = new Dictionary<int, List<double[]>>();
            Dictionary<int, List<double[]>> gPoints = new Dictionary<int, List<double[]>>();
            Dictionary<int, List<double[]>> nodesStress = new Dictionary<int, List<double[]>>();
            Dictionary<int, List<double[]>> nodesStrain = new Dictionary<int, List<double[]>>();
            //for (int i = 1; i < allStepsSolutions.Count; i++)
            for (int i = 0; i <= allStepsSolutions.Count - 1; i++)
            {
                elementsAssembly.UpdateDisplacements(allStepsSolutions[i]);
                if(i == allStepsSolutions.Count - 1)
                {
                    stress = elementsAssembly.GetElementsStresses(allStepsSolutions[i]);
                    strain = elementsAssembly.GetElementsStains(allStepsSolutions[i]);
                    gPoints = elementsAssembly.GetElementsGaussPoints(allStepsSolutions[i]);
                    nodesStress = elementsAssembly.GetElementsNodesStresses(allStepsSolutions[i]);
                    nodesStrain = elementsAssembly.GetElementsNodesStains(allStepsSolutions[i]);
                    ExportToFile.ExportMatlabFinalGeometry(elementsAssembly, BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[i], elementsAssembly.BoundedDOFsVector));
                }
            }
            ShowToGUI.PlotFinalGeometry(elementsAssembly);
            double[] fullDynamicSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[200], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[210], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[220], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[230], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[240], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol6 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[250], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol7 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[260], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol8 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[270], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol9 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[280], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol10 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[285], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol11 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[290], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol12 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[291], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol13 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[292], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol14 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[293], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol15 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[294], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol16 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[295], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol17 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[296], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol18 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[297], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol19 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[298], elementsAssembly.BoundedDOFsVector);
            double[] fullDynamicSol20 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[299], elementsAssembly.BoundedDOFsVector);
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
            VectorOperations.PrintVectorToFile(fullDynamicSol11, @"C:\Users\Public\Documents\Results11.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol12, @"C:\Users\Public\Documents\Results12.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol13, @"C:\Users\Public\Documents\Results13.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol14, @"C:\Users\Public\Documents\Results14.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol15, @"C:\Users\Public\Documents\Results15.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol16, @"C:\Users\Public\Documents\Results16.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol17, @"C:\Users\Public\Documents\Results17.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol18, @"C:\Users\Public\Documents\Results18.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol19, @"C:\Users\Public\Documents\Results19.dat");
            VectorOperations.PrintVectorToFile(fullDynamicSol20, @"C:\Users\Public\Documents\Results20.dat");
            VectorOperations.PrintDictionaryofListsofVectorsToFile(stress, @"C:\Users\Public\Documents\Stress.dat");
            VectorOperations.PrintDictionaryofListsofVectorsToFile(strain, @"C:\Users\Public\Documents\Strain.dat");
            VectorOperations.PrintDictionaryofListsofVectorsToFile(gPoints, @"C:\Users\Public\Documents\GaussPoints.dat");
            VectorOperations.PrintDictionaryofListsofVectorsToFile(nodesStress, @"C:\Users\Public\Documents\StressNodes.dat");
            VectorOperations.PrintDictionaryofListsofVectorsToFile(nodesStrain, @"C:\Users\Public\Documents\StrainNodes.dat");
            //newSolver.PrintExplicitSolution();
            Results finalResults = new Results() { DynamicSolution = newSolver.explicitSolution, TimeSteps = newSolver.TimeAtEachStep, SelectedDOF = 1, SelectedInterval = 1, SolutionType = "Dynamic" };
            return finalResults;
        }
    }
}
