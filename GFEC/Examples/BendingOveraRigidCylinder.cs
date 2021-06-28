using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class BendingOveraRigidCylinder
    {
        public static ISolver structuralSolution;
        static int[] structuralBoundaryConditions;
        const double thickness = 0.3;
        const double thickness2 = 0.4;
        const double gap = 0.001;
        private const int circles = 3;
        private const int steps = 36;
        private const double Radius = 0.50;
        private const double radius = 0.35;
        private const double RD = (Radius - radius) / (circles - 1);
        private const double angle = (Math.PI / 180) * (360 / steps);
        const int nodesInX = 25;
        const int nodesInY = 2;
        const int nodesInΖ = 2;
        const int nodesNumber = 316;
        const int elementsNumber = 96;
        const int contactElements = 20;
        const double xInterv = 0.25;
        const double yInterv = 0.25;

        //External loads
        const double externalStructuralLoad = -200.0;

        static List<int> loadedStructuralDOFs;
        static double[] externalForcesStructuralVector;

        const double YoungMod = 1.0 * 1e8;
        static double[] center = { -3.0, -Radius - gap, 0.0 };

        const double poissonRatio = 0.0;
        const double density = 8000.0;
        const double area = 1.0;
        const double contactArea = thickness * xInterv / 4.0;



        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            for(int i = 1; i <= 12; i++)
            {
                boundedDofs.Add(3* i -2);
                boundedDofs.Add(3* i - 1);
                boundedDofs.Add(3* i);
            }
            for (int i = 101; i <= nodesNumber; i++)
            {
                boundedDofs.Add(3 * i - 2);
                boundedDofs.Add(3 * i - 1);
                boundedDofs.Add(3 * i);
            }
            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            loadedStructuralDOFs.Add(296);
            loadedStructuralDOFs.Add(299);
            externalForcesStructuralVector = new double[nodesNumber * 3];
        }

        private static Dictionary<int, INode> CreateNodes()
        {

            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            int l;
            l = 1;
            //First beam
            for (int i = 0; i < nodesInX; i++)
            {
                for (int j = 0; j < nodesInY; j++)
                {
                    for(int k = 0; k< nodesInΖ; k++)
                    {
                        nodes[l] = new Node( -i * xInterv, j * yInterv, k * thickness - thickness/2.0);
                        l += 1;
                    }
                }
            }
            //Rigid cylinder
            for (int i = circles; i > 0; i--)
            {
                for (int j = 0; j < steps; j++)
                {
                    nodes[l] = new Node(center[0] + (radius + RD * (i - 1)) * Math.Cos(j * (angle)), center[1] + ((radius + RD * (i - 1)) * Math.Sin(j * (angle))), -thickness2/2.0);
                    l += 1;
                }
            }
            for (int i = circles; i > 0; i--)
            {
                for (int j = 0; j < steps; j++)
                {
                    nodes[l] = new Node(center[0] + (radius + RD * (i - 1)) * Math.Cos(j * (angle)), center[1] + ((radius + RD * (i - 1)) * Math.Sin(j * (angle))), thickness2/2.0);
                    l += 1;
                }
            }
            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int k = 1;
            int n1 = 5;
            for (int i = 1; i <= nodesInX - 1; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, n1 }, { 2, n1 - 4 }, { 3, n1 - 2 },
                    { 4, n1 + 2 },{ 5, n1 + 1 },{ 6, n1 - 3 },{ 7, n1 - 1 },{ 8, n1 + 3 } };
                k += 1;
                n1 += 4;
            }
            for (int i = 1; i < steps; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i + 100 }, { 2, i + 101 },
                        { 3, i + 1 + steps + 100 }, { 4, i + steps + 100 },
                        { 5, i + 208}, { 6, i + 209},
                        { 7, i + 1 + steps + 208}, { 8, i + steps + 208}};
                k += 1;
            }
            connectivity[k] = new Dictionary<int, int>() { { 1, steps + 100 }, { 2, 101},
                    { 3, steps + 101  }, { 4, 2 * steps + 100 },
                    { 5, steps + 208 }, { 6, 209},
                    { 7, steps + 209 }, { 8, 2 * steps + 208 } };// close circle
            k += 1;
            for (int i = 1; i < steps; i++)
            {
                connectivity[k] = new Dictionary<int, int>() { { 1, i + steps + 100 }, { 2, i + 101 + steps },
                        { 3, i + 101 + 2 * steps}, { 4, i + 2 * steps + 100 },
            { 5, i + steps + 208 }, { 6, i + 209 + steps },
                        { 7, i + 209 + 2 * steps}, { 8, i + 2 * steps + 208 }};
                k += 1;
            }
            connectivity[k] = new Dictionary<int, int>() { { 1, 2 * steps + 100 }, { 2, steps + 101},
                    { 3, 2 * steps + 101 }, { 4, 3 * steps + 100 },
            { 5, 2 * steps + 208 }, { 6, steps + 209},
                    { 7, 2 * steps + 209 }, { 8, 3 * steps + 208 } };// close circle

            k += 1;
            //Contact elements
            connectivity[k] = new Dictionary<int, int>() { { 1, 118 }, { 2, 226 }, { 3, 225 }, {4, 117 }, {5, 57} };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 118 }, { 2, 226 }, { 3, 225 }, { 4, 117 }, { 5, 58 } };
            k += 1; 
            connectivity[k] = new Dictionary<int, int>() { { 1, 117 }, { 2, 225 }, { 3, 224 }, { 4, 116 }, { 5, 57 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 117 }, { 2, 225 }, { 3, 224 }, { 4, 116 }, { 5, 58 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 116 }, { 2, 224 }, { 3, 223 }, { 4, 115 }, { 5, 57 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 116 }, { 2, 224 }, { 3, 223 }, { 4, 115 }, { 5, 58 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 115 }, { 2, 223 }, { 3, 222 }, { 4, 114 }, { 5, 57 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 115 }, { 2, 223 }, { 3, 222 }, { 4, 114 }, { 5, 58 } };
            k += 1;

            connectivity[k] = new Dictionary<int, int>() { { 1, 114 }, { 2, 222 }, { 3, 221 }, { 4, 113 }, { 5, 53 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 114 }, { 2, 222 }, { 3, 221 }, { 4, 113 }, { 5, 54 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 113 }, { 2, 221 }, { 3, 220 }, { 4, 112 }, { 5, 53 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 113 }, { 2, 221 }, { 3, 220 }, { 4, 112 }, { 5, 54 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 112 }, { 2, 220 }, { 3, 219 }, { 4, 111 }, { 5, 53 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 112 }, { 2, 220 }, { 3, 219 }, { 4, 111 }, { 5, 54 } };
            k += 1;

            connectivity[k] = new Dictionary<int, int>() { { 1, 110 }, { 2, 218 }, { 3, 217 }, { 4, 109 }, { 5, 49 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 110 }, { 2, 218 }, { 3, 217 }, { 4, 109 }, { 5, 50 } };
            k += 1;

            connectivity[k] = new Dictionary<int, int>() { { 1, 108 }, { 2, 216 }, { 3, 215 }, { 4, 107 }, { 5, 45 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 108 }, { 2, 216 }, { 3, 215 }, { 4, 107 }, { 5, 46 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 107 }, { 2, 215 }, { 3, 214 }, { 4, 106 }, { 5, 45 } };
            k += 1;
            connectivity[k] = new Dictionary<int, int>() { { 1, 107 }, { 2, 215 }, { 3, 214 }, { 4, 106 }, { 5, 46 } };
            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= nodesNumber; i++)
            {
                nodeFAT[i] = new bool[] { true, true, true, false, false, false };
            }
            return nodeFAT;
        }
        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;

            double A = area;
            double CA = contactArea;

            string type = "Hex8";
            string type2 = "ContactNtS3D";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= 24; i++)
            {
                elementProperties[i] = new ElementProperties(E, poissonRatio, A, thickness, density, type);
            }
            for (int i = 25; i <= elementsNumber; i++)
            {
                elementProperties[i] = new ElementProperties(100.0 * E, poissonRatio, A, thickness, density, type);
            }
            for (int i = elementsNumber + 1; i <= elementsNumber + contactElements; i++)
            {
                elementProperties[i] = new ElementProperties(E, CA, type2);
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
            //ShowToGUI.PlotInitialGeometry(elementsAssembly);
            structuralSolution.LinearScheme = new LUFactorization();
            structuralSolution.NonLinearScheme.Tolerance = 1e-5;
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 10;

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
            //Dictionary<int, List<double[]>> nodalStress = new Dictionary<int, List<double[]>>();
            //Dictionary<int, List<double[]>> nodalStrain = new Dictionary<int, List<double[]>>();
            for (int i = 1; i <= allStepsSolutions.Count; i++)
            {
                string name = "NodalCoordinates" + i.ToString() + ".dat";
                //ExportToFile.ExportUpdatedNodalCoordinates(elementsAssembly, BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions.Single(m => m.Key == i).Value, elementsAssembly.BoundedDOFsVector), name);
                gPointsStress = elementsAssembly.GetElementsStresses(allStepsSolutions[i]);
                gPointsStrain = elementsAssembly.GetElementsStains(allStepsSolutions[i]);
                gPoints = elementsAssembly.GetElementsGaussPoints(allStepsSolutions[i]);
                //nodalStress = elementsAssembly.GetElementsNodesStresses(allStepsSolutions[i]);
                //nodalStrain = elementsAssembly.GetElementsNodesStains(allStepsSolutions[i]);
                string name1 = "GPointsStress" + i.ToString() + ".dat";
                string name2 = "GPointsStrain" + i.ToString() + ".dat";
                string name3 = "GPointsCoordinates" + i.ToString() + ".dat";
                //string name4 = "NodalStress" + i.ToString() + ".dat";
                //string name5 = "NodalStrain" + i.ToString() + ".dat";

                VectorOperations.PrintDictionaryofListsofVectorsToFile(gPointsStress, @"C:\Users\Public\Documents\" + name1);
                VectorOperations.PrintDictionaryofListsofVectorsToFile(gPointsStrain, @"C:\Users\Public\Documents\" + name2);
                VectorOperations.PrintDictionaryofListsofVectorsToFile(gPoints, @"C:\Users\Public\Documents\" + name3);
                //VectorOperations.PrintDictionaryofListsofVectorsToFile(nodalStress, @"C:\Users\Public\Documents\" + name4);
                //VectorOperations.PrintDictionaryofListsofVectorsToFile(nodalStrain, @"C:\Users\Public\Documents\" + name5);
            }
            elementsAssembly.UpdateDisplacements(solvector3);
            //ShowToGUI.PlotFinalGeometry(elementsAssembly);
            double[] fullSolVector3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solvector3, elementsAssembly.BoundedDOFsVector);
            //Dictionary<int, INode> finalNodes = Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullSolVector3);
            //double[] xFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item1;
            //double[] yFinalNodalCoor = Assembly.NodalCoordinatesToVectors(finalNodes).Item2;
            Dictionary<int, double[]> allStepsFullSolutions = new Dictionary<int, double[]>();
            Dictionary<int, Dictionary<int, double[]>> allStepsContactForces = new Dictionary<int, Dictionary<int, double[]>>();
            Dictionary<int, double[]> elementsInternalContactForcesVector;
            for (int i = 1; i <= allStepsSolutions.Count; i++)
            {
                elementsInternalContactForcesVector = new Dictionary<int, double[]>();
                elementsAssembly.UpdateDisplacements(allStepsSolutions[i]);
                for (int j = 1; j <= contactElements; j++)
                {
                    elementsInternalContactForcesVector[elementsNumber + j] = elementsAssembly.ElementsAssembly[elementsNumber + j].CreateInternalGlobalForcesVector();
                }
                allStepsContactForces[i] = elementsInternalContactForcesVector;
                string name = "ContactForces" + i.ToString() + ".dat";
                double[] Vector = new double[60];
                int count = 0;
                for(int j = 0; j < 20; j++)
                {
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j + 1).Value[12];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j + 1).Value[13];
                    count += 1;
                    Vector[count] = allStepsContactForces.Single(m => m.Key == i).Value.Single(n => n.Key == elementsNumber + j + 1).Value[14];
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
    }
}
