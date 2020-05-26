using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;

namespace GFEC
{
    public static class CoupledThermalStructuralCNTsInAngle3
    {
        private const double Lenght = 8.0;
        private const double Overlap = 0.12;
        private const int totalNodes = 426;
        private const int totalContactElements = 10;//20;//8;
        private const int totalElements = 280;
        private const int nodesInXCoor = 81;
        private const int nodesInXCoor2 = 61;
        private const int nodesInYCoor = 3;
        private const double scaleFactor = 1.0;
        private const double xIntervals = 0.1;
        private const double xIntervals2 = 0.2;
        private const double yIntervals = 0.1;
        private const double ContactLength = Overlap * Lenght;
        private const double offset = Lenght - ContactLength;//- 0.10 * xIntervals;//8.1;//9.3//-0.65 * xIntervals;
        //private const double gap = 1.139;
        private const double gap = 0.381;
        //private const double angle = Math.PI / 2.2;
        private const double angle = Math.PI * 0.48485;
        private const int loadDof1 = (nodesInXCoor * (nodesInYCoor - 1) + nodesInXCoor - totalContactElements + 1) * 2;
        private const int LoadDof2 = loadDof1 + (totalContactElements - 1) * 2;
        private const int ThermalDof1 = 2;
        private const int ThermalDof2 = nodesInXCoor * (nodesInYCoor - 1) + 2;
        public static ISolver structuralSolution;
        //public static ISolver thermalSolution;

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
        const double Pressure = -20.00 * 1000.0;
        const double T0 = 100.0;
        const double cond = 3300 * 1.0e-6;
        const double thickness = 0.1;
        const double externalStructuralLoad = Pressure * thickness * xIntervals;
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
        const double density = 1400.0;
        const double area = thickness * xIntervals;
        const double solidThermalCond = cond;
        //const double roughness = (1.0 / 2.81);
        const double roughness = 2.81e-6;
        const double contactCond = cond;
        const double yieldStrength = 60.0 * 1.0e6;
        //----------------------------------------------------------------------
        //const double YoungMod = 1.0e-6;
        //const double density = 8000.0;
        //const double area = 0.01;
        //const double thickness = 0.1;
        //const double solidThermalCond = 3300 * 1.0e-9;
        //const double roughness = (1.0 / 2.81);
        //const double contactCond = 3300 * 1.0e-9;
        //const double yieldStrength = 60.0 * 1e-9;


        //private static void CreateStructuralBoundaryConditions()
        //{
        //    List<int> boundedDofs = new List<int>();
        //    for (int i = 0; i < nodesInYCoor; i++)
        //    {
        //        boundedDofs.Add(i * 2 * nodesInXCoor + 1); //upper beam left side support
        //        boundedDofs.Add(i * 2 * nodesInXCoor + 2); //upper beam left side support
        //    }
        //    for (int i = 1; i <= totalContactElements; i++)
        //    {
        //        boundedDofs.Add(nodesInXCoor * nodesInYCoor * 2 + 2 * i); //lower beam lower side support
        //    }
        //    for (int i = 1; i <= nodesInYCoor; i++)
        //    {
        //        boundedDofs.Add(nodesInYCoor * nodesInXCoor * 2 + nodesInXCoor2 * 2 * i - 1); //lower beam right side support
        //        boundedDofs.Add(nodesInYCoor * nodesInXCoor * 2 + nodesInXCoor2 * 2 * i); //lower beam right side support

        //    }
        //    //boundedDofs.Add(2);//Additional restraints to prevent Rigid Body Motion
        //    //boundedDofs.Add(totalNodes * 2);
        //    structuralBoundaryConditions = boundedDofs.ToArray<int>();
        //}
        private static void CreateStructuralBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            for (int i = 0; i < nodesInYCoor; i++) //upper beam left side support
            {
                boundedDofs.Add(i * nodesInXCoor * 2 + 1);
                boundedDofs.Add(i * nodesInXCoor * 2 + 2);
            }
            for (int i = nodesInXCoor * nodesInYCoor + 1; i < totalNodes; i++)
            {
                boundedDofs.Add(i * 2 + 0);
                boundedDofs.Add(i * 2 + 1); //lower beam support for all nodes
            }

            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateThermalBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();

            for (int i = 1; i <= nodesInYCoor; i++)
            {
                boundedDofs.Add(nodesInYCoor * nodesInXCoor + nodesInXCoor2 * i); //lower beam right side support
            }
            for (int i = 0; i < nodesInYCoor; i++)
            {
                boundedDofs.Add(i * nodesInXCoor + 1); //upper beam left side support
            }
            thermalBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            for (int i = 0; i < (3/2) * totalContactElements; i++)
            {
                loadedStructuralDOFs.Add(nodesInXCoor * nodesInYCoor * 2 - 2 * i);
            }
            externalForcesStructuralVector = new double[totalNodes * 2];
        }
        //private static void CreateStructuralLoadVector()
        //{
        //    loadedStructuralDOFs = new List<int>();
        //    for (int i = 0; i < nodesInXCoor - 1; i++)
        //    {
        //        loadedStructuralDOFs.Add(nodesInXCoor * nodesInYCoor * 2 - 2 * i);
        //    }
        //    externalForcesStructuralVector = new double[totalNodes * 2];
        //}

        private static void CreateThermalLoadVector()
        {
            loadedThermalDOFs = new List<int>();
            for (int i = 0; i < nodesInYCoor; i++)
            {
                loadedThermalDOFs.Add(nodesInXCoor * i + 2);
            }
            externalHeatLoafVector = new double[totalNodes];
            //thermalBoundaryConditions2 = loadedThermalDOFs.ToArray<int>();
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
                for (int i = 0; i < nodesInXCoor2; i++)
                {
                    nodes[k] = new Node(i * xIntervals2 * scaleFactor + offset, j * yIntervals * scaleFactor - ((nodesInYCoor - 1) * yIntervals + gap));
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

            for (int j = 0; j <= nodesInYCoor - 2; j++)
            {
                for (int i = 1; i <= nodesInXCoor2 - 1; i++)
                {
                    connectivity[k] = new Dictionary<int, int>() { { 1, i + nodesInYCoor * nodesInXCoor + j * nodesInXCoor2 }, { 2, i + 1 + nodesInYCoor * nodesInXCoor + j * nodesInXCoor2 }, { 3, i + nodesInXCoor2 + 1 + nodesInYCoor * nodesInXCoor + j * nodesInXCoor2 }, { 4, i + nodesInXCoor2 + nodesInYCoor * nodesInXCoor + j * nodesInXCoor2 } };
                    k += 1;
                }

            }

            //Contact elements
            k = 1;
            for (int i = 1; i <= totalContactElements / 2; i++)
            {
                //int lowerMiddleNode =  (nodesInXCoor + nodesInXCoor2) * nodesInYCoor - nodesInXCoor2 + i + 1;
                int lowerLeftNode = (nodesInXCoor + nodesInXCoor2) * nodesInYCoor - nodesInXCoor2 + i;
                int lowerRightNode = lowerLeftNode + 1;
                int upperNode = nodesInXCoor - totalContactElements + k;
                connectivity[totalElements + k] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
                connectivity[totalElements + k + 1] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode + 1 } };
                k = k + 2;

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

            for (int i = totalElements + 1; i <= totalElements + totalContactElements; i++)
            {
                if (i == totalElements + totalContactElements)
                {
                    elementProperties[i] = new ElementProperties(E, A / 2, type2);

                }
                else
                {
                    elementProperties[i] = new ElementProperties(E, A, type2);
                }
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
            string type3 = "Quad4Th2";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= (nodesInXCoor - 1) * (nodesInYCoor - 1); i++)
            {
                elementProperties[i] = new ElementProperties();
                elementProperties[i].ElementType = type;
                elementProperties[i].ThermalConductivity = thermalCond;
            }
            for (int i = 1 + (nodesInXCoor - 1) * (nodesInYCoor - 1); i <= totalElements; i++)
            {
                elementProperties[i] = new ElementProperties();
                elementProperties[i].ElementType = type3;
                elementProperties[i].ThermalConductivity = thermalCond;
                elementProperties[i].A = xIntervals2;
                elementProperties[i].B = yIntervals;
            }
            for (int i = totalElements + 1; i <= totalElements + totalContactElements; i++)
            {
                elementProperties[i] = new ElementProperties();
                elementProperties[i].ElementType = type2;
                elementProperties[i].ThermalConductivity = thermalCond;
                if (i == totalElements + totalContactElements)
                {
                    elementProperties[i].SectionArea = A / 2;

                }
                else
                {
                    elementProperties[i].SectionArea = A;
                }
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
            //string pathForContour1 = @"C:\Users\Public\Documents\Total\1";
            //string pathForContour2 = @"C:\Users\Public\Documents\Total\2";
            //ExportToFile.CreateContourDataForMatlab(Xvec1Initial, Yvec1Initial, Ζvec1Initial, nodesInYCoor, nodesInXCoor, pathForContour1);
            //ExportToFile.CreateContourDataForMatlab(Xvec2Initial, Yvec2Initial, Ζvec2Initial, nodesInYCoor, nodesInXCoor, pathForContour2);
            //ISolver structuralSolution = new StaticSolver();
            structuralSolution.LinearScheme = new PCGSolver();
            //structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
            //structuralSolution.NonLinearScheme.Tolerance = 1e-5;
            structuralSolution.NonLinearScheme.Tolerance = 0.00001;
            //structuralSolution.NonLinearScheme.MaxIterations = 100;
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 20;

            double[] externalForces3 = externalForcesStructuralVector;
            foreach (var dof in loadedStructuralDOFs)
            {
                if (dof == loadDof1 | dof == LoadDof2)
                {
                    externalForces3[dof - 1] = externalStructuralLoad / 2;

                }
                else
                {
                    externalForces3[dof - 1] = externalStructuralLoad;

                }
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
            Dictionary<int, Dictionary<int, double>> allStepsContactDXs = new Dictionary<int, Dictionary<int, double>>();
            Dictionary<int, Dictionary<int, double>> allStepsContactDXs1 = new Dictionary<int, Dictionary<int, double>>();
            Dictionary<int, Dictionary<int, double>> allStepsContactDXs2 = new Dictionary<int, Dictionary<int, double>>();
            Dictionary<int, Dictionary<int, INode>> allstepsfinalNodes = new Dictionary<int, Dictionary<int, INode>>();
            //Dictionary<int, Dictionary<int, INode>> PreviusStatefinalNodes = new Dictionary<int, Dictionary<int, INode>>();
            Dictionary<int, double[]> elementsInternalContactForcesVector;
            Dictionary<int, double> Deltax;
            Dictionary<int, double> Deltax1;
            Dictionary<int, double> Deltax2;
            //double[] SolutionSum = new double[2 * totalNodes];
            double[] NodalX = new double[totalNodes];
            double[] NodalY = new double[totalNodes];
            int upperNode = new int();
            int lowerLeftNode = new int();
            int lowerRightNode = new int();
            for (int i = 1; i <= allStepsSolutions.Count; i++)
            {
                elementsInternalContactForcesVector = new Dictionary<int, double[]>();
                Deltax = new Dictionary<int, double>();
                Deltax1 = new Dictionary<int, double>();
                Deltax2 = new Dictionary<int, double>();
                //allstepsfinalNodes = new Dictionary<int, Dictionary<int, INode>>();
                double[] fullSolVectorSteps = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[i], elementsAssembly.BoundedDOFsVector);
                //SolutionSum = VectorOperations.VectorVectorAddition(SolutionSum, fullSolVectorSteps);
                allstepsfinalNodes[i] = Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullSolVectorSteps);
                elementsAssembly.UpdateDisplacements(allStepsSolutions[i]);
                NodalX = Assembly.NodalCoordinatesToVectors(allstepsfinalNodes[i]).Item1;
                NodalY = Assembly.NodalCoordinatesToVectors(allstepsfinalNodes[i]).Item2;
                for (int j = totalElements + 1; j <= totalElements + totalContactElements; j++)
                {
                    elementsInternalContactForcesVector[j] = elementsAssembly.ElementsAssembly[j].CreateInternalGlobalForcesVector();
                }
                upperNode = nodesInXCoor - totalContactElements + 1;
                lowerLeftNode = (nodesInXCoor + nodesInXCoor2) * nodesInYCoor - nodesInXCoor2 + 1;
                lowerRightNode = lowerLeftNode + 1;
                for (int j = totalElements + 1; j <= totalElements + totalContactElements; j += 2)
                {
                    Deltax[j] = Math.Pow(Math.Pow(NodalX[lowerRightNode - 1] - NodalX[lowerLeftNode - 1], 2) + Math.Pow(NodalY[lowerRightNode - 1] - NodalY[lowerLeftNode - 1], 2), 0.5);
                    Deltax2[j] = Math.Pow(Math.Pow(NodalX[lowerRightNode - 1] - NodalX[upperNode - 1], 2) + Math.Pow(NodalY[lowerRightNode - 1] - NodalY[upperNode - 1], 2), 0.5);
                    Deltax1[j] = Math.Pow(Math.Pow(NodalX[upperNode - 1] - NodalX[lowerLeftNode - 1], 2) + Math.Pow(NodalY[upperNode - 1] - NodalY[lowerLeftNode - 1], 2), 0.5);
                    upperNode = upperNode + 1;
                    Deltax[j + 1] = Math.Pow(Math.Pow(NodalX[lowerRightNode - 1] - NodalX[lowerLeftNode - 1], 2) + Math.Pow(NodalY[lowerRightNode - 1] - NodalY[lowerLeftNode - 1], 2), 0.5);
                    Deltax2[j + 1] = Math.Pow(Math.Pow(NodalX[lowerRightNode - 1] - NodalX[upperNode - 1], 2) + Math.Pow(NodalY[lowerRightNode - 1] - NodalY[upperNode - 1], 2), 0.5);
                    Deltax1[j + 1] = Math.Pow(Math.Pow(NodalX[upperNode - 1] - NodalX[lowerLeftNode - 1], 2) + Math.Pow(NodalY[upperNode - 1] - NodalY[lowerLeftNode - 1], 2), 0.5);
                    lowerLeftNode = lowerRightNode;
                    lowerRightNode = lowerRightNode + 1;
                    upperNode = upperNode + 1;
                }
                allStepsContactForces[i] = elementsInternalContactForcesVector;
                allStepsContactDXs[i] = Deltax;
                allStepsContactDXs1[i] = Deltax1;
                allStepsContactDXs2[i] = Deltax2;

            }
            List<double[]> structuralSolutions = new List<double[]>();

            #endregion
            #region Thermal
            List<double[]> thermalSolutions = new List<double[]>();
            //List<double[]> thermalSolutions2 = new List<double[]>();
            List<Dictionary<int, double>> contactContactivityForEachStep = new List<Dictionary<int, double>>();
            List<double[]> contactContactivityForEachStep2 = new List<double[]>();
            for (int k = 1; k <= allStepsSolutions.Count; k++)
            {
                IAssembly elementsAssembly2 = CreateThermalAssembly();

                for (int j = totalElements + 1; j <= totalElements + totalContactElements; j++)
                {
                    double[] contactForce = allStepsContactForces[k][j];
                    double DX = allStepsContactDXs[k][j];
                    double DX1 = allStepsContactDXs1[k][j];
                    double DX2 = allStepsContactDXs2[k][j];
                    //double ContactForceVar = Math.Pow(Math.Pow(VectorOperations.VectorNorm2(new double[] { contactForce[4], contactForce[5] }),2),0.5);
                    elementsAssembly2.ElementsProperties[j].ContactForceValue = VectorOperations.VectorNorm2(new double[] { contactForce[4], contactForce[5] });
                    elementsAssembly2.ElementsProperties[j].Dx = DX;
                    elementsAssembly2.ElementsProperties[j].Dx1 = DX1;
                    elementsAssembly2.ElementsProperties[j].Dx2 = DX2;

                }

                elementsAssembly2.CreateElementsAssembly();
                elementsAssembly2.ActivateBoundaryConditions = true;
                double[,] globalStiffnessMatrix2 = elementsAssembly2.CreateTotalStiffnessMatrix();

                ISolver thermalSolution = new StaticSolver();
                thermalSolution.LinearScheme = new LUFactorization();
                thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                thermalSolution.NonLinearScheme.Tolerance = 1e-7;
                thermalSolution.NonLinearScheme.MaxIterations = 100;
                thermalSolution.ActivateNonLinearSolver = true;
                thermalSolution.NonLinearScheme.numberOfLoadSteps = 10;

                thermalSolution.AssemblyData = elementsAssembly2;
                double[] externalHeatFlux = externalHeatLoafVector;
                //externalHeatFlux[0] =  externalHeatLoad;
                //externalHeatFlux[15] = externalHeatLoad;
                //externalHeatFlux[30] = externalHeatLoad;
                //externalHeatFlux[45] = externalHeatLoad;
                //externalHeatFlux[60] = externalHeatLoad;

                foreach (var dof in loadedThermalDOFs)
                {
                    if ((dof == ThermalDof1 | dof == ThermalDof2))
                    {
                        externalHeatFlux[dof - 1] = -T0 * (cond / (6 * xIntervals * yIntervals)) * ((Math.Pow(xIntervals, 2) - 2 * Math.Pow(yIntervals, 2)) - (Math.Pow(xIntervals, 2) + Math.Pow(yIntervals, 2)));
                    }
                    else
                    {
                        externalHeatFlux[dof - 1] = -2 * T0 * (cond / (6 * xIntervals * yIntervals)) * ((Math.Pow(xIntervals, 2) - 2 * Math.Pow(yIntervals, 2)) - (Math.Pow(xIntervals, 2) + Math.Pow(yIntervals, 2)));
                    }
                }
                //for (int i = 61; i <= 75; i++)
                //{
                //    externalHeatFlux[61] = externalHeatLoad;
                //}
                double[] reducedExternalHeatFlux = BoundaryConditionsImposition.ReducedVector(externalHeatFlux, thermalSolution.AssemblyData.BoundedDOFsVector);
                thermalSolution.Solve(reducedExternalHeatFlux);
                Dictionary<int, double[]> IntHeatFlux = thermalSolution.GetInternalForces();
                //double [] HeatSolve= new double[totalNodes];
                //double[] value;
                ////double[] Values = new double[totalNodes];
                //for (int i = 1; i <= 10; i++)
                //{
                //    value = IntHeatFlux[i];
                //    thermalSolutions2.Add(value);
                //    //Values[i - 1] = Convert.ToDouble(value);
                //}
                //HeatSolve = VectorOperations.VectorVectorAddition(Values , externalHeatFlux);
                double[] tempSol = thermalSolution.GetSolution();
                elementsAssembly2.UpdateDisplacements(tempSol);
                //double[] contactContactivityVector = new double[totalContactElements];
                //Dictionary<int, double[]> allStepsSolutions2 = thermalSolution.GetAllStepsSolutions();
                thermalSolutions.Add(tempSol);
                //thermalSolutions2.Add(HeatSolve);
                Dictionary<int, double> contactContactivity = AssemblyHelpMethods.RetrieveContactContactivity(thermalSolution.AssemblyData);
                //    contactContactivityForEachStep.Add(contactContactivity);
                //    for (int i = 0; i <= totalContactElements - 1; i++)
                //    {
                //        contactContactivityVector[i] = contactContactivity[i + totalElements + 1];
                //    }
                //    contactContactivityForEachStep2.Add(contactContactivityVector);
            }
            int[] thermalBoundCond = thermalBoundaryConditions;
            //int[] thermalBoundCond2 = thermalBoundaryConditions2;
            double[] fullStructuralSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[4], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[8], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[12], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[16], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[20], elementsAssembly.BoundedDOFsVector);
            double[] fullThermalSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[3], thermalBoundCond);
            double[] fullThermalSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[7], thermalBoundCond);
            double[] fullThermalSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[11], thermalBoundCond);
            double[] fullThermalSol4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[15], thermalBoundCond);
            double[] fullThermalSol5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[19], thermalBoundCond);
            //double[] fullThermalSol6 = thermalSolutions2[19];
            //double[] fullThermalSol7 = thermalSolutions2[39];
            //double[] fullThermalSol8 = thermalSolutions2[59];
            //double[] fullThermalSol9 = thermalSolutions2[79];
            //double[] fullThermalSol10 = thermalSolutions2[99];
            //double[] FullcontactContactivity1 = contactContactivityForEachStep2[1];
            //double[] FullcontactContactivity2 = contactContactivityForEachStep2[3];
            //double[] FullcontactContactivity3 = contactContactivityForEachStep2[5];
            //double[] FullcontactContactivity4 = contactContactivityForEachStep2[7];
            //double[] FullcontactContactivity5 = contactContactivityForEachStep2[9];
            double[] FullStructuralSolution_1_Χ = new double[totalNodes];
            double[] FullStructuralSolution_1_Υ = new double[totalNodes];
            double[] FullStructuralSolution_2_Χ = new double[totalNodes];
            double[] FullStructuralSolution_2_Υ = new double[totalNodes];
            double[] FullStructuralSolution_3_Χ = new double[totalNodes];
            double[] FullStructuralSolution_3_Υ = new double[totalNodes];
            double[] FullStructuralSolution_4_Χ = new double[totalNodes];
            double[] FullStructuralSolution_4_Υ = new double[totalNodes];
            double[] FullStructuralSolution_5_Χ = new double[totalNodes];
            double[] FullStructuralSolution_5_Υ = new double[totalNodes];
            int count = 1;
            for (int i = 1; i <= 2 * totalNodes - 1; i += 2)
            {
                FullStructuralSolution_1_Χ[count - 1] = fullStructuralSol1[i - 1];
                FullStructuralSolution_1_Υ[count - 1] = fullStructuralSol1[i];
                count += 1;
            }
            count = 1;
            for (int i = 1; i <= 2 * totalNodes - 1; i += 2)
            {
                FullStructuralSolution_2_Χ[count - 1] = fullStructuralSol2[i - 1];
                FullStructuralSolution_2_Υ[count - 1] = fullStructuralSol2[i];
                count += 1;
            }
            count = 1;
            for (int i = 1; i <= 2 * totalNodes - 1; i += 2)
            {
                FullStructuralSolution_3_Χ[count - 1] = fullStructuralSol3[i - 1];
                FullStructuralSolution_3_Υ[count - 1] = fullStructuralSol3[i];
                count += 1;
            }
            count = 1;
            for (int i = 1; i <= 2 * totalNodes - 1; i += 2)
            {
                FullStructuralSolution_4_Χ[count - 1] = fullStructuralSol4[i - 1];
                FullStructuralSolution_4_Υ[count - 1] = fullStructuralSol4[i];
                count += 1;
            }
            count = 1;
            for (int i = 1; i <= 2 * totalNodes - 1; i += 2)
            {
                FullStructuralSolution_5_Χ[count - 1] = fullStructuralSol5[i - 1];
                FullStructuralSolution_5_Υ[count - 1] = fullStructuralSol5[i];
                count += 1;
            }
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol1), FullStructuralSolution_1_Χ, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.00Mpa_Res1_ux.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol2), FullStructuralSolution_2_Χ, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.00Mpa_Res2_ux.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol3), FullStructuralSolution_3_Χ, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.00Mpa_Res3_ux.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol4), FullStructuralSolution_4_Χ, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.00Mpa_Res4_ux.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol5), FullStructuralSolution_5_Χ, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.00Mpa_Res5_ux.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol1), FullStructuralSolution_1_Υ, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.00Mpa_Res1_uy.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol2), FullStructuralSolution_2_Υ, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.00Mpa_Res2_uy.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol3), FullStructuralSolution_3_Υ, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.00Mpa_Res3_uy.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol4), FullStructuralSolution_4_Υ, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.00Mpa_Res4_uy.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol5), FullStructuralSolution_5_Υ, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.00Mpa_Res5_uy.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol1), fullThermalSol1, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.000Mpa_Res1_temp.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol2), fullThermalSol2, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.000Mpa_Res2_temp.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol3), fullThermalSol3, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.000Mpa_Res3_temp.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol4), fullThermalSol4, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.000Mpa_Res4_temp.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol5), fullThermalSol5, @"C:\Users\Public\Documents\ThermalStructuralCNTsInAngle3_0.12L_20.000Mpa_Res5_temp.dat");
            //ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol1), fullThermalSol6, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.90L_0.100Mpa_Res6.dat");
            //ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol2), fullThermalSol7, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.90L_0.100Mpa_Res7.dat");
            //ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol3), fullThermalSol8, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.90L_0.100Mpa_Res8.dat");
            //ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol4), fullThermalSol9, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.90L_0.100Mpa_Res9.dat");
            //ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol5), fullThermalSol10, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.90L_0.100Mpa_Res10.dat");
            //VectorOperations.PrintVectorToFile(fullThermalSol6, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.30L_0.100Mpa_Res6.dat");
            //VectorOperations.PrintVectorToFile(fullThermalSol7, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.30L_0.100Mpa_Res7.dat");
            //VectorOperations.PrintVectorToFile(fullThermalSol8, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.30L_0.100Mpa_Res8.dat");
            //VectorOperations.PrintVectorToFile(fullThermalSol9, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.30L_0.100Mpa_Res9.dat");
            //VectorOperations.PrintVectorToFile(fullThermalSol10, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.30L_0.100Mpa_Res10.dat");

            //VectorOperations.PrintVectorToFile(FullcontactContactivity1, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.30L_0.100Mpa_CC1.dat");
            //VectorOperations.PrintVectorToFile(FullcontactContactivity2, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.30L_0.100Mpa_CC2.dat");
            //VectorOperations.PrintVectorToFile(FullcontactContactivity3, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.30L_0.100Mpa_CC3.dat");
            //VectorOperations.PrintVectorToFile(FullcontactContactivity4, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.30L_0.100Mpa_CC4.dat");
            //VectorOperations.PrintVectorToFile(FullcontactContactivity5, @"C:\Users\Public\Documents\CouplThermStructCNTs_0.30L_0.100Mpa_CC5.dat");
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

            //string path = @"C:\Users\Public\Documents\Total\1final";
            //string path2 = @"C:\Users\Public\Documents\Total\2final";
            //ExportToFile.CreateContourDataForMatlab(Xvec1Final, Yvec1Final, Ζvec1Final, nodesInYCoor, nodesInXCoor, path);
            //ExportToFile.CreateContourDataForMatlab(Xvec2Final, Yvec2Final, Ζvec2Final, nodesInYCoor, nodesInXCoor, path2);

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

        //public static void RunDynamicExample()
        //{
        //    IAssembly elementsAssembly = CreateAssembly();
        //    elementsAssembly.CreateElementsAssembly();
        //    elementsAssembly.ActivateBoundaryConditions = true;

        //    InitialConditions initialValues = new InitialConditions();
        //    initialValues.InitialAccelerationVector = new double[6];
        //    initialValues.InitialDisplacementVector = new double[6];
        //    //initialValues.InitialDisplacementVector[7] = -0.02146;
        //    initialValues.InitialVelocityVector = new double[6];
        //    initialValues.InitialTime = 0.0;

        //    ExplicitSolver newSolver = new ExplicitSolver(1.0, 10000);
        //    newSolver.Assembler = elementsAssembly;

        //    newSolver.InitialValues = initialValues;
        //    newSolver.ExternalForcesVector = new double[] { 0, 0, 0, 0, -50000, -50000 };
        //    newSolver.LinearSolver = new CholeskyFactorization();
        //    newSolver.ActivateNonLinearSolution = true;
        //    newSolver.SolveNewmark();
        //    newSolver.PrintExplicitSolution();//
        //}

    }
}

