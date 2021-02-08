using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;

namespace GFEC
{
    public static class CoupledThermalStructuralCNTsInAngle5
    {
        private const int totalNodes = 486;
        private const int totalContactElements = 30;//20;//8;
        private const int totalElements = 320;
        private const int nodesInXCoor = 81;
        private const int nodesInYCoor = 3;
        private const double scaleFactor = 1.0;
        private const double xIntervals = 0.1;
        private const double yIntervals = 0.1;
        private const double offset = 4.55;//8.1;//9.3// tested: 7.0 - 0.05, 7.0 -0.45, 7.0 -0.25, 7.0 -0.0
        private const double gap = 2.77; //tested: 1.14, 2.75, 2.10, 0.75
        public static ISolver structuralSolution;
        public static ISolver thermalSolution;
        private const double angle = Math.PI / 2.57; //tested: 2.2, 2.57, 2.40, 2.12
        private static int loadStepsNumber = 40;
        //private const double angle = Math.PI * 0.48485;

        //--------------------------------------------
        //private const int totalNodes = 150;
        //private const int totalContactElements = 8;
        //private const int totalElements = 112;
        //private const int nodesInXCoor = 15;
        //private const int nodesInYCoor = 5;
        //private const double scaleFactor = 1.0;
        //private const double xIntervals = 0.1;
        //private const double yIntervals = 0.1;
        //private const double offset = 0.7;
        //private const double gap = 0.01;



        //Boundary conditions
        //Model1
        //static readonly int[] structuralBoundaryConditions = new int[] { 1, 2, 31, 32, 61, 62, 91, 92, 121, 122, 152, 179, 180, 209, 210, 239, 240, 269, 270, 299, 300 };
        //static readonly int[] thermalBoundaryConditions = new int[] { 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90 };

        //Model2
        static int[] structuralBoundaryConditions; // = new int[] { 1, 203, 505, 707, 909, 1012, 1014, 1016, 1018, 1020, 1022, 1024, 1026, 1211, 1413, 1615, 1817, 2019 };
        static int[] thermalBoundaryConditions; //= new int[] { 606, 707, 808, 909, 1010 };




        //External loads
        const double externalStructuralLoad = -10000.0 * 2;
        const double externalHeatLoad = 2500.0 * 1e-6;
        //-----------------------------------------------------------------------------------
        //const double externalStructuralLoad = -5 * 100000000.0 * 1e-18 * 0.3;
        //const double externalHeatLoad = 2500.0 * 1e-9;






        static List<int> loadedStructuralDOFs; // = new List<int>(new int[] { 995, 997, 999, 1001, 1003, 1005, 1007, 1009 });
        static double[] externalForcesStructuralVector; // = new double[2020];



        static double[] externalHeatLoafVector; // = new double[1010];
        //static readonly List<int> loadedThermalDOFs = new List<int>(new int[] { 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75 });
        static List<int> loadedThermalDOFs; // = new List<int>(new int[] { 0, 101, 202, 303, 404 }); //zero indexed

        //CNT values
        //const double YoungMod = 1.0e12;
        //const double density = 8000.0;
        //const double area = 0.01;
        //const double thickness = 0.1;
        //const double solidThermalCond = 6600;
        //const double roughness = (1000000.0 / 2.81);
        //const double contactCond = 6600;
        //const double yieldStrength = 60000000000.0;

        //CNT values scaled
        const double YoungMod = 1.0 * 1e9;
        const double density = 8000.0;
        const double area = 0.01;
        const double thickness = 0.1;
        const double solidThermalCond = 3300 * 1.0e-6;
        const double roughness = 2.81 * 1.0e-6;
        const double contactCond = 3300 * 1.0e-6;
        const double yieldStrength = 60.0 * 1e6;

        //----------------------------------------------------------------------
        //const double YoungMod = 1.0e-6;
        //const double density = 8000.0;
        //const double area = 0.01;
        //const double thickness = 0.1;
        //const double solidThermalCond = 3300 * 1.0e-9;
        //const double roughness = (1.0 / 2.81);
        //const double contactCond = 3300 * 1.0e-9;
        //const double yieldStrength = 60.0 * 1e-9;


        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            //for (int i = 0; i < nodesInYCoor; i++)
            //{
            //    boundedDofs.Add(i * 2 * nodesInXCoor + 1); //upper beam left side support
            //}

            //for (int i = 0; i < nodesInYCoor; i++)
            //{
            //    boundedDofs.Add(i * 2 * nodesInXCoor + 2 * nodesInXCoor - 1); //upper beam right side support
            //}

            for (int i = 0; i < nodesInYCoor; i++) //upper beam left side support
            {
                boundedDofs.Add(i * nodesInXCoor * 2 + 1);
                boundedDofs.Add(i * nodesInXCoor * 2 + 2);
            }

            //for (int i = 1; i <= totalContactElements; i++)
            //{
            //    boundedDofs.Add(nodesInXCoor * nodesInYCoor * 2 + 2 * i); //lower beam lower side support
            //}
            //for (int i = 0; i < nodesInYCoor; i++)
            //{
            //    boundedDofs.Add(nodesInYCoor * nodesInXCoor * 2 + nodesInXCoor * 2 * (i+1) - 1); //lower beam right side support
            //}

            //for (int i = 0; i < totalNodes; i++)
            //{
            //    boundedDofs.Add(i * 2 + 1); //support for all nodes at X direction
            //}

            for (int i = totalNodes / 2 + 1; i < totalNodes; i++)
            {
                boundedDofs.Add(i * 2 + 0);
                boundedDofs.Add(i * 2 + 1); //lower beam support for all nodes
            }

            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateThermalBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();

            for (int i = 0; i < nodesInYCoor; i++)
            {
                boundedDofs.Add(nodesInYCoor * nodesInXCoor + nodesInXCoor * (i + 1)); //lower beam right side support
            }
            thermalBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            for (int i = 0; i < totalContactElements; i++)
            {
                loadedStructuralDOFs.Add(nodesInXCoor * nodesInYCoor * 2 - 2 * i);
            }
            //for (int i = 1; i <= nodesInYCoor; i++)
            //{
            //    loadedStructuralDOFs.Add((nodesInXCoor * (nodesInYCoor - 1) + i) * 2); //load at upper surface of upper beam
            //}
            externalForcesStructuralVector = new double[totalNodes * 2];
        }

        private static void CreateThermalLoadVector()
        {
            loadedThermalDOFs = new List<int>();
            for (int i = 0; i < nodesInYCoor; i++)
            {
                loadedThermalDOFs.Add(nodesInXCoor * i + 1);
            }
            externalHeatLoafVector = new double[totalNodes];
        }

        private static Dictionary<int, INode> CreateNodes()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            //Upper cantilever
            int k;
            k = 1;
            for (int j = 0; j < nodesInYCoor; j++)
            {
                for (int i = 0; i < nodesInXCoor; i++)
                {
                    nodes[k] = new Node(i * xIntervals * scaleFactor * Math.Sin(angle) + j * yIntervals * scaleFactor * Math.Cos(angle), -Math.Cos(angle) * i * xIntervals * scaleFactor + Math.Sin(angle) * j * yIntervals * scaleFactor);//j * yIntervals * scaleFactor - (j * yIntervals * scaleFactor * Math.Cos(angle)));
                    k += 1;
                }
            }

            //Lower cantilever
            for (int j = 0; j < nodesInYCoor; j++)
            {
                for (int i = 0; i < nodesInXCoor; i++)
                {
                    nodes[k] = new Node(i * xIntervals * scaleFactor + offset, j * yIntervals * scaleFactor - ((nodesInYCoor - 1) * yIntervals + gap));
                    k += 1;
                }
            }
            return nodes;
        }

        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity()
        {

            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            int k = 1;
            for (int j = 0; j <= nodesInYCoor - 2; j++)
            {
                for (int i = 1; i <= nodesInXCoor - 1; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + j * nodesInXCoor }, { 2, i + 1 + j * nodesInXCoor }, { 3, i + 1 + nodesInXCoor + j * nodesInXCoor }, { 4, i + nodesInXCoor + j * nodesInXCoor } };
                    k += 1;
                }

            }

            for (int j = 0 + nodesInYCoor; j <= nodesInYCoor - 2 + nodesInYCoor; j++)
            {
                for (int i = 1; i <= nodesInXCoor - 1; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + j * nodesInXCoor }, { 2, i + 1 + j * nodesInXCoor }, { 3, i + 1 + nodesInXCoor + j * nodesInXCoor }, { 4, i + nodesInXCoor + j * nodesInXCoor } };
                    k += 1;
                }

            }

            //Contact elements
            for (int i = 1; i <= totalContactElements - 1; i++)
            {
                int lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i + 1;
                int lowerLeftNode = lowerMiddleNode - 1;
                int lowerRightNode = lowerMiddleNode + 1;
                int upperNode = nodesInXCoor - totalContactElements + i + 1;
                connectivity[totalElements + i] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
            }
            for (int i = 1; i <= totalContactElements - 2; i++)
            {
                int lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i + 1;
                int lowerLeftNode = lowerMiddleNode - 1;
                int lowerRightNode = lowerMiddleNode + 1;
                int upperNode = nodesInXCoor - totalContactElements + i + 2;
                connectivity[totalElements + i + totalContactElements -1] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
            }
            for (int i = 1; i <= totalContactElements - 3; i++)
            {
                int lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i + 1;
                int lowerLeftNode = lowerMiddleNode - 1;
                int lowerRightNode = lowerMiddleNode + 1;
                int upperNode = nodesInXCoor - totalContactElements + i + 3;
                connectivity[totalElements + i + totalContactElements - 1 + totalContactElements - 2] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
            }

            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= totalNodes; i++)
            {
                nodeFAT[i] = new bool[] { true, true, false, false, false, false };
            }
            return nodeFAT;
        }

        private static Dictionary<int, bool[]> CreateThermalNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= totalNodes; i++)
            {
                nodeFAT[i] = new bool[] { true, false, false, false, false, false };
            }
            return nodeFAT;
        }

        private static Dictionary<int, IElementProperties> CreateElementProperties()
        {
            double E = YoungMod;
            double A = area;
            string type = "Quad4";
            string type2 = "ContactNtS2D";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= totalElements; i++)
            {
                elementProperties[i] = new ElementProperties(E, A, type);
            }

            for (int i = 1; i <= totalElements; i++)
            {
                elementProperties[i].Density = density;
                elementProperties[i].Thickness = thickness;
            }

            for (int i = totalElements + 1; i <= totalElements + totalContactElements - 1 + totalContactElements -2 + totalContactElements - 3; i++)
            {
                elementProperties[i] = new ElementProperties(E, A, type2);
                elementProperties[i].Density = density;
                elementProperties[i].Thickness = thickness;
            }
            return elementProperties;
        }

        private static Dictionary<int, IElementProperties> CreateThermalElementProperties()
        {
            double thermalCond = solidThermalCond;
            double A = area;
            string type = "Quad4Th";
            string type2 = "ContactNtS2DTh";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= totalElements; i++)
            {
                elementProperties[i] = new ElementProperties();
                elementProperties[i].ElementType = type;
                elementProperties[i].ThermalConductivity = thermalCond;
            }
            for (int i = totalElements + 1; i <= totalElements + totalContactElements - 1 + totalContactElements - 2+ totalContactElements - 3; i++)
            {
                elementProperties[i] = new ElementProperties();
                elementProperties[i].ElementType = type2;
                elementProperties[i].ThermalConductivity = thermalCond;
                elementProperties[i].SectionArea = A;
                elementProperties[i].SurfaceRoughness = roughness;
                elementProperties[i].ContactThermalConductivity = contactCond;
                elementProperties[i].YieldStrength = yieldStrength;
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

        private static IAssembly CreateThermalAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes();
            assembly.ElementsConnectivity = CreateConnectivity();
            assembly.ElementsProperties = CreateThermalElementProperties();
            assembly.NodeFreedomAllocationList = CreateThermalNodeFAT();
            CreateThermalBoundaryConditions();
            CreateThermalLoadVector();
            assembly.BoundedDOFsVector = thermalBoundaryConditions;
            return assembly;
        }

        public static Results RunStaticExample()
        {
            #region Structural
            IAssembly elementsAssembly = CreateAssembly();
            elementsAssembly.CreateElementsAssembly();
            elementsAssembly.ActivateBoundaryConditions = true;
            double[,] globalStiffnessMatrix = elementsAssembly.CreateTotalStiffnessMatrix();

            //Gnuplot graphs
            ShowToGUI.PlotInitialGeometry(elementsAssembly);

            Dictionary<int, INode> initialNodes = elementsAssembly.Nodes;
            double[] initialXCoord = Assembly.NodalCoordinatesToVectors(initialNodes).Item1;
            double[] initialYCoord = Assembly.NodalCoordinatesToVectors(initialNodes).Item2;

            double[] Xvec1Initial = new double[totalNodes / 2];
            double[] Yvec1Initial = new double[totalNodes / 2];
            double[] Xvec2Initial = new double[totalNodes / 2];
            double[] Yvec2Initial = new double[totalNodes / 2];
            double[] Ζvec1Initial = Enumerable.Repeat(1.0, totalNodes / 2).ToArray();
            double[] Ζvec2Initial = Enumerable.Repeat(1.0, totalNodes / 2).ToArray();

            Array.Copy(initialXCoord, 0, Xvec1Initial, 0, totalNodes / 2);
            Array.Copy(initialYCoord, 0, Yvec1Initial, 0, totalNodes / 2);

            Array.Copy(initialXCoord, totalNodes / 2, Xvec2Initial, 0, totalNodes / 2);
            Array.Copy(initialYCoord, totalNodes / 2, Yvec2Initial, 0, totalNodes / 2);
            string pathForContour1 = @"C:\Users\Public\Documents\Total\1";
            string pathForContour2 = @"C:\Users\Public\Documents\Total\2";
            ExportToFile.CreateContourDataForMatlab(Xvec1Initial, Yvec1Initial, Ζvec1Initial, nodesInYCoor, nodesInXCoor, pathForContour1);
            ExportToFile.CreateContourDataForMatlab(Xvec2Initial, Yvec2Initial, Ζvec2Initial, nodesInYCoor, nodesInXCoor, pathForContour2);




            ///structuralSolution = new StaticSolver();
            structuralSolution.LinearScheme = new PCGSolver();
            //structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
            structuralSolution.NonLinearScheme.Tolerance = 1e-5;
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.numberOfLoadSteps = loadStepsNumber;

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

            Dictionary<int, Dictionary<int, double[]>> allStepsContactForces = new Dictionary<int, Dictionary<int, double[]>>();
            Dictionary<int, Dictionary<int, double>> allStepsProjectionPoints = new Dictionary<int, Dictionary<int, double>>();
            Dictionary<int, double[]> elementsInternalContactForcesVector;
            Dictionary<int, double> projectionPointForEachElement;
            for (int i = 1; i <= allStepsSolutions.Count; i++)
            {
                elementsInternalContactForcesVector = new Dictionary<int, double[]>();
                projectionPointForEachElement = new Dictionary<int, double>();
                elementsAssembly.UpdateDisplacements(allStepsSolutions[i]);
                for (int j = totalElements + 1; j <= totalElements + totalContactElements - 1; j++)
                {
                    elementsInternalContactForcesVector[j] = elementsAssembly.ElementsAssembly[j].CreateInternalGlobalForcesVector();
                    projectionPointForEachElement[j] = elementsAssembly.ElementsAssembly[j].ClosestPointProjection();
                }
                allStepsContactForces[i] = elementsInternalContactForcesVector;
                allStepsProjectionPoints[i] = projectionPointForEachElement;
            }



            List<double[]> structuralSolutions = new List<double[]>();

            ExportToFile.ExportMatlabInitialGeometry(elementsAssembly);
            #endregion


            #region Thermal
            Dictionary<int, double[]> thermalSolutions = new Dictionary<int, double[]>();
            Dictionary<int, Dictionary<int, double[]>> allStepsHeatFluxes = new Dictionary<int, Dictionary<int, double[]>>();
            List<Dictionary<int, double>> contactContactivityForEachStep = new List<Dictionary<int, double>>();
            for (int k = 1; k <= allStepsSolutions.Count; k++)
            {
                IAssembly elementsAssembly2 = CreateThermalAssembly();

                for (int j = totalElements + 1; j <= totalElements + totalContactElements - 1; j++)
                {
                    double[] contactForce = allStepsContactForces[k][j];
                    elementsAssembly2.ElementsProperties[j].ContactForceValue = -contactForce[5];
                    double projectionPoint = allStepsProjectionPoints[k][j];
                    elementsAssembly2.ElementsProperties[j].Dx1 = projectionPoint;
                }

                elementsAssembly2.CreateElementsAssembly();
                elementsAssembly2.ActivateBoundaryConditions = true;
                double[,] globalStiffnessMatrix2 = elementsAssembly2.CreateTotalStiffnessMatrix();

                //ISolver thermalSolution = new StaticSolver();
                thermalSolution.LinearScheme = new LUFactorization();
                //thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                thermalSolution.NonLinearScheme.Tolerance = 1e-7;
                thermalSolution.ActivateNonLinearSolver = true;
                thermalSolution.NonLinearScheme.numberOfLoadSteps = 20;

                thermalSolution.AssemblyData = elementsAssembly2;
                double[] externalHeatFlux = externalHeatLoafVector;
                //externalHeatFlux[0] =  externalHeatLoad;
                //externalHeatFlux[15] = externalHeatLoad;
                //externalHeatFlux[30] = externalHeatLoad;
                //externalHeatFlux[45] = externalHeatLoad;
                //externalHeatFlux[60] = externalHeatLoad;

                foreach (var dof in loadedThermalDOFs)
                {
                    externalHeatFlux[dof - 1] = externalHeatLoad;
                }
                //for (int i = 61; i <= 75; i++)
                //{
                //    externalHeatFlux[61] = externalHeatLoad;
                //}
                double[] reducedExternalHeatFlux = BoundaryConditionsImposition.ReducedVector(externalHeatFlux, thermalSolution.AssemblyData.BoundedDOFsVector);
                thermalSolution.Solve(reducedExternalHeatFlux);
                double[] tempSol = thermalSolution.GetSolution();
                thermalSolutions.Add(k, tempSol);

                double[] fullThermalSolutionVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(tempSol, elementsAssembly2.BoundedDOFsVector);
                elementsAssembly2.UpdateDisplacements(fullThermalSolutionVector);
                Dictionary<int, double[]> elementsInternalHeatFluxesVector = new Dictionary<int, double[]>();
                for (int j = totalElements + 1; j <= totalElements + totalContactElements - 1; j++)
                {
                    elementsInternalHeatFluxesVector[j] = elementsAssembly2.ElementsAssembly[j].CreateInternalGlobalForcesVector();
                }
                allStepsHeatFluxes[k] = elementsInternalHeatFluxesVector;

                Dictionary<int, double> contactContactivity = AssemblyHelpMethods.RetrieveContactContactivity(thermalSolution.AssemblyData);
                contactContactivityForEachStep.Add(contactContactivity);

            }

            ExportToFile.ExportGeometryDataWithTemperatures(structuralSolution, thermalSolutions, thermalBoundaryConditions);
            ExportToFile.ExportCondactivityForAllLoadSteps(contactContactivityForEachStep);
            ExportToFile.ExportContactForcesForAllLoadSteps(allStepsContactForces);
            ExportToFile.ExportHeatFluxesForAllLoadSteps(allStepsHeatFluxes);
            ExportToFile.ExportConvergenceResultsToFile(structuralSolution.NonLinearScheme.LoadStepConvergence);

            int[] thermalBoundCond = thermalBoundaryConditions;
            double[] fullStructuralSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[8], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[16], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[24], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[32], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[40], elementsAssembly.BoundedDOFsVector);
            double[] fullThermalSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[7], thermalBoundCond);
            double[] fullThermalSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[15], thermalBoundCond);
            double[] fullThermalSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[23], thermalBoundCond);
            double[] fullThermalSol4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[31], thermalBoundCond);
            double[] fullThermalSol5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[39], thermalBoundCond);
            double[] contactContactivityForLoadStep1 = contactContactivityForEachStep[7].Values.ToArray();
            double[] contactContactivityForLoadStep2 = contactContactivityForEachStep[15].Values.ToArray();
            double[] contactContactivityForLoadStep3 = contactContactivityForEachStep[23].Values.ToArray();
            double[] contactContactivityForLoadStep4 = contactContactivityForEachStep[31].Values.ToArray();
            double[] contactContactivityForLoadStep5 = contactContactivityForEachStep[39].Values.ToArray();
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol1), fullThermalSol1, @"C:\Users\Public\Documents\Results1.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol2), fullThermalSol2, @"C:\Users\Public\Documents\Results2.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol3), fullThermalSol3, @"C:\Users\Public\Documents\Results3.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol4), fullThermalSol4, @"C:\Users\Public\Documents\Results4.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol5), fullThermalSol5, @"C:\Users\Public\Documents\Results5.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep1, @"C:\Users\Public\Documents\contactivity1.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep2, @"C:\Users\Public\Documents\contactivity2.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep3, @"C:\Users\Public\Documents\contactivity3.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep4, @"C:\Users\Public\Documents\contactivity4.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep5, @"C:\Users\Public\Documents\contactivity5.dat");

            structuralSolutions.Add(fullStructuralSol1);
            structuralSolutions.Add(fullStructuralSol2);
            structuralSolutions.Add(fullStructuralSol3);
            structuralSolutions.Add(fullStructuralSol4);
            structuralSolutions.Add(fullStructuralSol5);



            double[] Xvec1Final = new double[totalNodes / 2];
            double[] Yvec1Final = new double[totalNodes / 2];
            double[] Xvec2Final = new double[totalNodes / 2];
            double[] Yvec2Final = new double[totalNodes / 2];
            double[] Ζvec1Final = new double[totalNodes / 2];
            double[] Ζvec2Final = new double[totalNodes / 2];

            Array.Copy(xFinalNodalCoor, 0, Xvec1Final, 0, totalNodes / 2);
            Array.Copy(yFinalNodalCoor, 0, Yvec1Final, 0, totalNodes / 2);
            Array.Copy(fullThermalSol4, 0, Ζvec1Final, 0, totalNodes / 2);
            Array.Copy(xFinalNodalCoor, totalNodes / 2, Xvec2Final, 0, totalNodes / 2);
            Array.Copy(yFinalNodalCoor, totalNodes / 2, Yvec2Final, 0, totalNodes / 2);
            Array.Copy(fullThermalSol4, totalNodes / 2, Ζvec2Final, 0, totalNodes / 2);

            List<HeatMapData> plots2 = new List<HeatMapData>();
            plots2.Add(new HeatMapData() { Xcoordinates = Xvec1Final, Ycoordinates = Yvec1Final, Temperatures = Ζvec1Final });
            plots2.Add(new HeatMapData() { Xcoordinates = Xvec2Final, Ycoordinates = Yvec2Final, Temperatures = Ζvec2Final });

            ShowToGUI.PlotHeatMap(plots2);

            string path = @"C:\Users\Public\Documents\Total\1final";
            string path2 = @"C:\Users\Public\Documents\Total\2final";
            ExportToFile.CreateContourDataForMatlab(Xvec1Final, Yvec1Final, Ζvec1Final, nodesInYCoor, nodesInXCoor, path);
            ExportToFile.CreateContourDataForMatlab(Xvec2Final, Yvec2Final, Ζvec2Final, nodesInYCoor, nodesInXCoor, path2);

            //ExportToFile.ExportGeometryDataWithTemperatures(finalNodes, fullTempSol);

            GnuPlot.Close();

            while (true)
            {
                if (File.Exists(AppContext.BaseDirectory + "gnuplot.png") && new FileInfo(AppContext.BaseDirectory + "gnuplot.png").Length > 0)
                {
                    break;

                }
                Thread.Sleep(100);
            }
            GnuPlot.KillProcess();
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

