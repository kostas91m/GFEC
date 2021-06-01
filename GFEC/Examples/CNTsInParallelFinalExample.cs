using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;

namespace GFEC
{
    public static class CNTsInParallelFinalExample
    {
        private const int offsetNodes = 73; /*73;*///65;/*55;*//*41;*//*29;*///17;//2
        private const int totalNodes = 648;
        private const int addedNodesrodElements = 2;
        private const int totalElements = 480;
        private const int nodesInXCoor = 81;
        private const int nodesInYCoor = 4;
        private const int totalContactElements = nodesInXCoor - offsetNodes + 1;//20;//8;
        private const double scaleFactor = 1.0;
        private const double xIntervals = 0.375;
        private const double yIntervals = 0.41;
        private const double offset = (offsetNodes-1)*xIntervals;//8.1;//9.3;
        private const double gap = 0.00000001;
        public static ISolver structuralSolution;

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
        private const int ThermalDof1 = 2;
        private const int ThermalDof2 = nodesInXCoor * (nodesInYCoor - 1) + 2;
        const double externalStructuralLoad = -2.85 * 81/totalContactElements;
        const double T0 = 100.0;
        const double solidThermalCond = 3300;
        static double externalHeatLoad = -2 * T0 * (solidThermalCond / (6 * xIntervals * yIntervals)) * ((Math.Pow(xIntervals, 2) - 2 * Math.Pow(yIntervals, 2)) - (Math.Pow(xIntervals, 2) + Math.Pow(yIntervals, 2)));        //const double externalStructuralLoad = -5 * 100000000.0 * 1e-18 * 0.3;
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
        const double YoungMod = 1.45e6;
        const double density = 8000.0;

        const double thickness = 0.38;
        const double area = thickness * yIntervals;
        //const double roughness = 2.81 * 1.0e-6;
        const double roughness = 0.0075;
        const double contactCond = 3300;
        const double yieldStrength = 60.0 * 1e3;

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
            for (int i = 0; i < nodesInYCoor; i++)
            {
                boundedDofs.Add(i * 2 * nodesInXCoor + 1); //upper beam left side support
            }
            for (int i = 0; i < nodesInYCoor; i++)
            {
                boundedDofs.Add(i * 2 * nodesInXCoor + 2 * nodesInXCoor - 1); //upper beam right side support
            }

            for (int i = 1; i <= addedNodesrodElements; i++) //upper beam rod elements supports
            {
                boundedDofs.Add(totalNodes * 2 + 2 * i - 1);
                boundedDofs.Add(totalNodes * 2 + 2 * i);
            }
            for (int i = 1; i <= totalContactElements; i++)
            {
                boundedDofs.Add(nodesInXCoor * nodesInYCoor * 2 + 2 * i); //lower beam lower side support
            }
            for (int i = 0; i < nodesInYCoor; i++)
            {
                boundedDofs.Add(nodesInYCoor * nodesInXCoor * 2 + nodesInXCoor * 2 * (i + 1) - 1); //lower beam right side support
                boundedDofs.Add(nodesInYCoor * nodesInXCoor * 2 + nodesInXCoor * 2 * (i + 1)); //lower beam right side support
            }

            //for (int i = 0; i < totalNodes; i++)
            //{
            //    boundedDofs.Add(i * 2 + 1); //support for all nodes at X direction
            //}

            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateThermalBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();

            for (int i = 0; i < nodesInYCoor; i++)
            {
                boundedDofs.Add(nodesInYCoor * nodesInXCoor + nodesInXCoor * (i + 1)); //lower beam right side support
            }

            for (int i = 0; i < nodesInYCoor; i++)
            {
                boundedDofs.Add(nodesInXCoor * i + 1); //upper beam left side support
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
            externalForcesStructuralVector = new double[(totalNodes + addedNodesrodElements) * 2];
            //externalForcesStructuralVector = new double[(totalNodes) * 2];

        }

        private static void CreateThermalLoadVector()
        {
            loadedThermalDOFs = new List<int>();
            for (int i = 0; i < nodesInYCoor; i++)
            {
                loadedThermalDOFs.Add(nodesInXCoor * i + 2);
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
                    nodes[k] = new Node(i * xIntervals * scaleFactor, j * yIntervals * scaleFactor);
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
            //Rod elements added nodes
            //for (int j = 0; j < addedNodesrodElements; j++)
            //{
            //    nodes[totalNodes + j + 1] = new Node(j * xIntervals * scaleFactor * (nodesInXCoor - 1), yIntervals * nodesInYCoor);
            //}
            nodes[totalNodes + 1] = new Node(offset, yIntervals * nodesInYCoor);
            nodes[totalNodes + 2] = new Node(xIntervals * scaleFactor * (nodesInXCoor - 1), yIntervals * nodesInYCoor);
            //nodes[totalNodes + addedNodesrodElements + 2] = new Node((totalContactElements - 1) * xIntervals, yIntervals * nodesInYCoor);

            return nodes;
        }
        private static Dictionary<int, INode> CreateNodes2()
        {
            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            //Upper cantilever
            int k;
            k = 1;
            for (int j = 0; j < nodesInYCoor; j++)
            {
                for (int i = 0; i < nodesInXCoor; i++)
                {
                    nodes[k] = new Node(i * xIntervals * scaleFactor, j * yIntervals * scaleFactor);
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
            for (int i = 1; i <= totalContactElements; i++)
            {
                int lowerNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i;
                int upperNode = nodesInXCoor - totalContactElements + i;
                connectivity[totalElements + i] = new Dictionary<int, int>() { { 1, lowerNode }, { 2, upperNode } };
            }
            ////Rod elements
            //connectivity[totalElements + totalContactElements + 1] = new Dictionary<int, int>() { { 1, nodesInXCoor * (nodesInYCoor - 1) + 1 }, { 2, totalNodes + 1 } };
            connectivity[totalElements + totalContactElements + 1] = new Dictionary<int, int>() { { 1, nodesInXCoor * (nodesInYCoor - 1) + offsetNodes }, { 2, totalNodes + 1 } };
            connectivity[totalElements + totalContactElements + 2] = new Dictionary<int, int>() { { 1, nodesInXCoor * nodesInYCoor }, { 2, totalNodes + 2 } };
            //connectivity[totalElements + totalContactElements + 4] = new Dictionary<int, int>() { { 1, nodesInXCoor * (nodesInYCoor - 1) + totalContactElements }, { 2, totalNodes + 4 } };
            return connectivity;
        }
        private static Dictionary<int, Dictionary<int, int>> CreateConnectivity2()
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
            for (int i = 1; i <= totalContactElements; i++)
            {
                int lowerNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i;
                int upperNode = nodesInXCoor - totalContactElements + i;
                connectivity[totalElements + i] = new Dictionary<int, int>() { { 1, lowerNode }, { 2, upperNode } };
            }
            return connectivity;
        }
        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= totalNodes + addedNodesrodElements; i++)
            //for (int i = 1; i <= totalNodes; i++)
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
            double contactElemArea = thickness * xIntervals;
            //double rodElementsYoungMod = 50 * 1e3;
            double rodElementsYoungMod = 4800;
            double rodElementsSectionArea = 0.01 / 4;
            string type = "Quad4";
            string type2 = "ContactNtN2D";
            string type3 = "Bar2D";
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
                if(i != totalElements + 1 && i!= totalElements + totalContactElements)
                {
                    elementProperties[i] = new ElementProperties(E, contactElemArea, type2);
                    elementProperties[i].Density = density;
                    elementProperties[i].Thickness = thickness;
                }
                else
                {
                    elementProperties[i] = new ElementProperties(E, contactElemArea/2, type2);
                    elementProperties[i].Density = density;
                    elementProperties[i].Thickness = thickness;
                }

            }
            //Rod elemnts properties
            for (int i = totalElements + totalContactElements + 1; i <= totalElements + totalContactElements + addedNodesrodElements; i++)
            {
                elementProperties[i] = new ElementProperties(rodElementsYoungMod/10, rodElementsSectionArea/10, type3);
            }
            return elementProperties;
        }
        private static Dictionary<int, IElementProperties> CreateThermalElementProperties()
        {
            double thermalCond = solidThermalCond;
            double A = thickness * xIntervals;
            string type = "Quad4Th";
            string type2 = "ContactNtN2DTh";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= totalElements; i++)
            {
                elementProperties[i] = new ElementProperties();
                elementProperties[i].ElementType = type;
                elementProperties[i].ThermalConductivity = thermalCond;
            }
            for (int i = totalElements + 1; i <= totalElements + totalContactElements; i++)
            {
                if (i != totalElements + 1 && i != totalElements + totalContactElements)
                {
                    elementProperties[i] = new ElementProperties();
                    elementProperties[i].ElementType = type2;
                    elementProperties[i].ThermalConductivity = thermalCond;
                    elementProperties[i].SectionArea = A;
                    elementProperties[i].SurfaceRoughness = roughness;
                    elementProperties[i].ContactThermalConductivity = contactCond;
                    elementProperties[i].YieldStrength = yieldStrength;
                }
                else
                {
                    elementProperties[i] = new ElementProperties();
                    elementProperties[i].ElementType = type2;
                    elementProperties[i].ThermalConductivity = thermalCond;
                    elementProperties[i].SectionArea = A/2;
                    elementProperties[i].SurfaceRoughness = roughness;
                    elementProperties[i].ContactThermalConductivity = contactCond;
                    elementProperties[i].YieldStrength = yieldStrength;
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

        private static IAssembly CreateThermalAssembly()
        {
            IAssembly assembly = new Assembly();
            assembly.Nodes = CreateNodes2();
            assembly.ElementsConnectivity = CreateConnectivity2();
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
            structuralSolution.LinearScheme = new LUFactorization();
            //structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
            structuralSolution.NonLinearScheme.Tolerance = 1e-5;
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 40;

            double[] externalForces3 = externalForcesStructuralVector;
            foreach (var dof in loadedStructuralDOFs)
            {
                if (dof == 2 * nodesInXCoor * nodesInYCoor || dof == 2 * (nodesInXCoor * nodesInYCoor - totalContactElements + 1))
                {
                    externalForces3[dof - 1] = externalStructuralLoad/2;

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
            Dictionary<int, double[]> elementsInternalContactForcesVector;
            for (int i = 1; i <= allStepsSolutions.Count; i++)
            {
                elementsInternalContactForcesVector = new Dictionary<int, double[]>();
                elementsAssembly.UpdateDisplacements(allStepsSolutions[i]);
                for (int j = totalElements + 1; j <= totalElements + totalContactElements; j++)
                {
                    elementsInternalContactForcesVector[j] = elementsAssembly.ElementsAssembly[j].CreateInternalGlobalForcesVector();
                }
                allStepsContactForces[i] = elementsInternalContactForcesVector;
            }



            List<double[]> structuralSolutions = new List<double[]>();

            #endregion


            #region Thermal
            Dictionary<int, double[]> thermalSolutions = new Dictionary<int, double[]>();
            Dictionary<int, Dictionary<int, double[]>> allStepsHeatFluxes = new Dictionary<int, Dictionary<int, double[]>>();
            List<Dictionary<int, double>> contactContactivityForEachStep = new List<Dictionary<int, double>>();
            for (int k = 1; k <= allStepsSolutions.Count; k++)
            {
                IAssembly elementsAssembly2 = CreateThermalAssembly();

                for (int j = totalElements + 1; j <= totalElements + totalContactElements; j++)
                {
                    double[] contactForce = allStepsContactForces[k][j];
                    elementsAssembly2.ElementsProperties[j].ContactForceValue = VectorOperations.VectorNorm2(new double[] { contactForce[2], contactForce[3] });
                }
                elementsAssembly2.CreateElementsAssembly();
                elementsAssembly2.ActivateBoundaryConditions = true;
                double[,] globalStiffnessMatrix2 = elementsAssembly2.CreateTotalStiffnessMatrix();

                ISolver thermalSolution = new StaticSolver();
                thermalSolution.LinearScheme = new LUFactorization();
                thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                thermalSolution.NonLinearScheme.Tolerance = 1e-9;
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
                        externalHeatFlux[dof - 1] = externalHeatLoad / 2;
                    }
                    else
                    {
                        externalHeatFlux[dof - 1] = externalHeatLoad;
                    }
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
            //ExportToFile.ExportGeometryDataWithTemperatures(structuralSolution, thermalSolutions, thermalBoundaryConditions);
            //ExportToFile.ExportCondactivityForAllLoadSteps(contactContactivityForEachStep);
            //ExportToFile.ExportContactForcesForAllLoadSteps(allStepsContactForces);
            //ExportToFile.ExportHeatFluxesForAllLoadSteps(allStepsHeatFluxes);

            int[] thermalBoundCond = thermalBoundaryConditions;
            double[] fullStructuralSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[1], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[2], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[3], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[4], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[5], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol6 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[6], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol7 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[7], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol8 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[8], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol9 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[9], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol10 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[10], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol11 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[11], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol12 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[12], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol13 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[13], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol14 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[14], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol15 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[15], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol16 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[16], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol17 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[17], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol18 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[18], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol19 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[19], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol20 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[20], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol21 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[21], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol22 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[22], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol23 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[23], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol24 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[24], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol25 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[25], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol26 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[26], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol27 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[27], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol28 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[28], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol29 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[29], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol30 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[30], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol31 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[31], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol32 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[32], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol33 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[33], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol34 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[34], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol35 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[35], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol36 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[36], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol37 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[37], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol38 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[38], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol39 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[39], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol40 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[40], elementsAssembly.BoundedDOFsVector);
            double[] fullThermalSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[1], thermalBoundCond);
            double[] fullThermalSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[2], thermalBoundCond);
            double[] fullThermalSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[3], thermalBoundCond);
            double[] fullThermalSol4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[4], thermalBoundCond);
            double[] fullThermalSol5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[5], thermalBoundCond);
            double[] fullThermalSol6 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[6], thermalBoundCond);
            double[] fullThermalSol7 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[7], thermalBoundCond);
            double[] fullThermalSol8 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[8], thermalBoundCond);
            double[] fullThermalSol9 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[9], thermalBoundCond);
            double[] fullThermalSol10 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[10], thermalBoundCond);
            double[] fullThermalSol11 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[11], thermalBoundCond);
            double[] fullThermalSol12 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[12], thermalBoundCond);
            double[] fullThermalSol13 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[13], thermalBoundCond);
            double[] fullThermalSol14 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[14], thermalBoundCond);
            double[] fullThermalSol15 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[15], thermalBoundCond);
            double[] fullThermalSol16 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[16], thermalBoundCond);
            double[] fullThermalSol17 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[17], thermalBoundCond);
            double[] fullThermalSol18 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[18], thermalBoundCond);
            double[] fullThermalSol19 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[19], thermalBoundCond);
            double[] fullThermalSol20 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[20], thermalBoundCond);
            double[] fullThermalSol21 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[21], thermalBoundCond);
            double[] fullThermalSol22 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[22], thermalBoundCond);
            double[] fullThermalSol23 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[23], thermalBoundCond);
            double[] fullThermalSol24 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[24], thermalBoundCond);
            double[] fullThermalSol25 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[25], thermalBoundCond);
            double[] fullThermalSol26 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[26], thermalBoundCond);
            double[] fullThermalSol27 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[27], thermalBoundCond);
            double[] fullThermalSol28 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[28], thermalBoundCond);
            double[] fullThermalSol29 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[29], thermalBoundCond);
            double[] fullThermalSol30 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[30], thermalBoundCond);
            double[] fullThermalSol31 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[31], thermalBoundCond);
            double[] fullThermalSol32 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[32], thermalBoundCond);
            double[] fullThermalSol33 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[33], thermalBoundCond);
            double[] fullThermalSol34 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[34], thermalBoundCond);
            double[] fullThermalSol35 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[35], thermalBoundCond);
            double[] fullThermalSol36 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[36], thermalBoundCond);
            double[] fullThermalSol37 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[37], thermalBoundCond);
            double[] fullThermalSol38 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[38], thermalBoundCond);
            double[] fullThermalSol39 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[39], thermalBoundCond);
            double[] fullThermalSol40 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[40], thermalBoundCond);
            //double[] fullThermalSolfinal1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[1], thermalBoundCond);
            //double[] fullThermalSolfinal2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[2], thermalBoundCond);
            //double[] fullThermalSolfinal3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[3], thermalBoundCond);
            //double[] fullThermalSolfinal4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[4], thermalBoundCond);
            //double[] fullThermalSolfinal5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[5], thermalBoundCond);
            //double[] fullThermalSolfinal6 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[6], thermalBoundCond);
            //double[] fullThermalSolfinal7 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[7], thermalBoundCond);
            //double[] fullThermalSolfinal8 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[8], thermalBoundCond);
            //double[] fullThermalSolfinal9 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[9], thermalBoundCond);
            //double[] fullThermalSolfinal10 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[10], thermalBoundCond);
            //double[] fullThermalSolfinal11 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[11], thermalBoundCond);
            //double[] fullThermalSolfinal12 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[12], thermalBoundCond);
            //double[] fullThermalSolfinal13 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[13], thermalBoundCond);
            //double[] fullThermalSolfinal14 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[14], thermalBoundCond);
            //double[] fullThermalSolfinal15 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[15], thermalBoundCond);
            //double[] fullThermalSolfinal16 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[16], thermalBoundCond);
            //double[] fullThermalSolfinal17 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[17], thermalBoundCond);
            //double[] fullThermalSolfinal18 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[18], thermalBoundCond);
            //double[] fullThermalSolfinal19 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[19], thermalBoundCond);
            //double[] fullThermalSolfinal20 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[20], thermalBoundCond);
            //double[] fullThermalSolfinal21 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[21], thermalBoundCond);
            //double[] fullThermalSolfinal22 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[22], thermalBoundCond);
            //double[] fullThermalSolfinal23 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[23], thermalBoundCond);
            //double[] fullThermalSolfinal24 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[24], thermalBoundCond);
            //double[] fullThermalSolfinal25 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[25], thermalBoundCond);
            //double[] fullThermalSolfinal26 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[26], thermalBoundCond);
            //double[] fullThermalSolfinal27 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[27], thermalBoundCond);
            //double[] fullThermalSolfinal28 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[28], thermalBoundCond);
            //double[] fullThermalSolfinal29 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[29], thermalBoundCond);
            //double[] fullThermalSolfinal30 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[30], thermalBoundCond);
            //double[] fullThermalSolfinal31 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[31], thermalBoundCond);
            //double[] fullThermalSolfinal32 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[32], thermalBoundCond);
            //double[] fullThermalSolfinal33 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[33], thermalBoundCond);
            //double[] fullThermalSolfinal34 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[34], thermalBoundCond);
            //double[] fullThermalSolfinal35 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[35], thermalBoundCond);
            //double[] fullThermalSolfinal36 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[36], thermalBoundCond);
            //double[] fullThermalSolfinal37 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[37], thermalBoundCond);
            //double[] fullThermalSolfinal38 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[38], thermalBoundCond);
            //double[] fullThermalSolfinal39 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[39], thermalBoundCond);
            //double[] fullThermalSolfinal40 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[40], thermalBoundCond);
            double[] fullThermalSolfinal1 = new double[fullThermalSol1.Length + 2];
            double[] fullThermalSolfinal2 = new double[fullThermalSol2.Length + 2];
            double[] fullThermalSolfinal3 = new double[fullThermalSol3.Length + 2];
            double[] fullThermalSolfinal4 = new double[fullThermalSol4.Length + 2];
            double[] fullThermalSolfinal5 = new double[fullThermalSol5.Length + 2];
            double[] fullThermalSolfinal6 = new double[fullThermalSol1.Length + 2];
            double[] fullThermalSolfinal7 = new double[fullThermalSol2.Length + 2];
            double[] fullThermalSolfinal8 = new double[fullThermalSol3.Length + 2];
            double[] fullThermalSolfinal9 = new double[fullThermalSol4.Length + 2];
            double[] fullThermalSolfinal10 = new double[fullThermalSol5.Length + 2];
            double[] fullThermalSolfinal11 = new double[fullThermalSol1.Length + 2];
            double[] fullThermalSolfinal12 = new double[fullThermalSol2.Length + 2];
            double[] fullThermalSolfinal13 = new double[fullThermalSol3.Length + 2];
            double[] fullThermalSolfinal14 = new double[fullThermalSol4.Length + 2];
            double[] fullThermalSolfinal15 = new double[fullThermalSol5.Length + 2];
            double[] fullThermalSolfinal16 = new double[fullThermalSol1.Length + 2];
            double[] fullThermalSolfinal17 = new double[fullThermalSol2.Length + 2];
            double[] fullThermalSolfinal18 = new double[fullThermalSol3.Length + 2];
            double[] fullThermalSolfinal19 = new double[fullThermalSol4.Length + 2];
            double[] fullThermalSolfinal20 = new double[fullThermalSol5.Length + 2];
            double[] fullThermalSolfinal21 = new double[fullThermalSol1.Length + 2];
            double[] fullThermalSolfinal22 = new double[fullThermalSol2.Length + 2];
            double[] fullThermalSolfinal23 = new double[fullThermalSol3.Length + 2];
            double[] fullThermalSolfinal24 = new double[fullThermalSol4.Length + 2];
            double[] fullThermalSolfinal25 = new double[fullThermalSol5.Length + 2];
            double[] fullThermalSolfinal26 = new double[fullThermalSol1.Length + 2];
            double[] fullThermalSolfinal27 = new double[fullThermalSol2.Length + 2];
            double[] fullThermalSolfinal28 = new double[fullThermalSol3.Length + 2];
            double[] fullThermalSolfinal29 = new double[fullThermalSol4.Length + 2];
            double[] fullThermalSolfinal30 = new double[fullThermalSol5.Length + 2];
            double[] fullThermalSolfinal31 = new double[fullThermalSol1.Length + 2];
            double[] fullThermalSolfinal32 = new double[fullThermalSol2.Length + 2];
            double[] fullThermalSolfinal33 = new double[fullThermalSol3.Length + 2];
            double[] fullThermalSolfinal34 = new double[fullThermalSol4.Length + 2];
            double[] fullThermalSolfinal35 = new double[fullThermalSol5.Length + 2];
            double[] fullThermalSolfinal36 = new double[fullThermalSol1.Length + 2];
            double[] fullThermalSolfinal37 = new double[fullThermalSol2.Length + 2];
            double[] fullThermalSolfinal38 = new double[fullThermalSol3.Length + 2];
            double[] fullThermalSolfinal39 = new double[fullThermalSol4.Length + 2];
            double[] fullThermalSolfinal40 = new double[fullThermalSol5.Length + 2];
            for (int runs = 0; runs < fullThermalSol1.Length; runs++)
            {
                fullThermalSolfinal1[runs] = fullThermalSol1[runs];
                fullThermalSolfinal2[runs] = fullThermalSol2[runs];
                fullThermalSolfinal3[runs] = fullThermalSol3[runs];
                fullThermalSolfinal4[runs] = fullThermalSol4[runs];
                fullThermalSolfinal5[runs] = fullThermalSol5[runs];
                fullThermalSolfinal6[runs] = fullThermalSol6[runs];
                fullThermalSolfinal7[runs] = fullThermalSol7[runs];
                fullThermalSolfinal8[runs] = fullThermalSol8[runs];
                fullThermalSolfinal9[runs] = fullThermalSol9[runs];
                fullThermalSolfinal10[runs] = fullThermalSol10[runs];
                fullThermalSolfinal11[runs] = fullThermalSol11[runs];
                fullThermalSolfinal12[runs] = fullThermalSol12[runs];
                fullThermalSolfinal13[runs] = fullThermalSol13[runs];
                fullThermalSolfinal14[runs] = fullThermalSol14[runs];
                fullThermalSolfinal15[runs] = fullThermalSol15[runs];
                fullThermalSolfinal16[runs] = fullThermalSol16[runs];
                fullThermalSolfinal17[runs] = fullThermalSol17[runs];
                fullThermalSolfinal18[runs] = fullThermalSol18[runs];
                fullThermalSolfinal19[runs] = fullThermalSol19[runs];
                fullThermalSolfinal20[runs] = fullThermalSol20[runs];
                fullThermalSolfinal21[runs] = fullThermalSol21[runs];
                fullThermalSolfinal22[runs] = fullThermalSol22[runs];
                fullThermalSolfinal23[runs] = fullThermalSol23[runs];
                fullThermalSolfinal24[runs] = fullThermalSol24[runs];
                fullThermalSolfinal25[runs] = fullThermalSol25[runs];
                fullThermalSolfinal26[runs] = fullThermalSol26[runs];
                fullThermalSolfinal27[runs] = fullThermalSol27[runs];
                fullThermalSolfinal28[runs] = fullThermalSol28[runs];
                fullThermalSolfinal29[runs] = fullThermalSol29[runs];
                fullThermalSolfinal30[runs] = fullThermalSol30[runs];
                fullThermalSolfinal31[runs] = fullThermalSol31[runs];
                fullThermalSolfinal32[runs] = fullThermalSol32[runs];
                fullThermalSolfinal33[runs] = fullThermalSol33[runs];
                fullThermalSolfinal34[runs] = fullThermalSol34[runs];
                fullThermalSolfinal35[runs] = fullThermalSol35[runs];
                fullThermalSolfinal36[runs] = fullThermalSol36[runs];
                fullThermalSolfinal37[runs] = fullThermalSol37[runs];
                fullThermalSolfinal38[runs] = fullThermalSol38[runs];
                fullThermalSolfinal39[runs] = fullThermalSol39[runs];
                fullThermalSolfinal40[runs] = fullThermalSol40[runs];
            }
            double[] contactContactivityForLoadStep1 = contactContactivityForEachStep[0].Values.ToArray();
            double[] contactContactivityForLoadStep2 = contactContactivityForEachStep[1].Values.ToArray();
            double[] contactContactivityForLoadStep3 = contactContactivityForEachStep[2].Values.ToArray();
            double[] contactContactivityForLoadStep4 = contactContactivityForEachStep[3].Values.ToArray();
            double[] contactContactivityForLoadStep5 = contactContactivityForEachStep[4].Values.ToArray();
            double[] contactContactivityForLoadStep6 = contactContactivityForEachStep[5].Values.ToArray();
            double[] contactContactivityForLoadStep7 = contactContactivityForEachStep[6].Values.ToArray();
            double[] contactContactivityForLoadStep8 = contactContactivityForEachStep[7].Values.ToArray();
            double[] contactContactivityForLoadStep9 = contactContactivityForEachStep[8].Values.ToArray();
            double[] contactContactivityForLoadStep10 = contactContactivityForEachStep[9].Values.ToArray();
            double[] contactContactivityForLoadStep11 = contactContactivityForEachStep[10].Values.ToArray();
            double[] contactContactivityForLoadStep12 = contactContactivityForEachStep[11].Values.ToArray();
            double[] contactContactivityForLoadStep13 = contactContactivityForEachStep[12].Values.ToArray();
            double[] contactContactivityForLoadStep14 = contactContactivityForEachStep[13].Values.ToArray();
            double[] contactContactivityForLoadStep15 = contactContactivityForEachStep[14].Values.ToArray();
            double[] contactContactivityForLoadStep16 = contactContactivityForEachStep[15].Values.ToArray();
            double[] contactContactivityForLoadStep17 = contactContactivityForEachStep[16].Values.ToArray();
            double[] contactContactivityForLoadStep18 = contactContactivityForEachStep[17].Values.ToArray();
            double[] contactContactivityForLoadStep19 = contactContactivityForEachStep[18].Values.ToArray();
            double[] contactContactivityForLoadStep20 = contactContactivityForEachStep[19].Values.ToArray();
            double[] contactContactivityForLoadStep21 = contactContactivityForEachStep[20].Values.ToArray();
            double[] contactContactivityForLoadStep22 = contactContactivityForEachStep[21].Values.ToArray();
            double[] contactContactivityForLoadStep23 = contactContactivityForEachStep[22].Values.ToArray();
            double[] contactContactivityForLoadStep24 = contactContactivityForEachStep[23].Values.ToArray();
            double[] contactContactivityForLoadStep25 = contactContactivityForEachStep[24].Values.ToArray();
            double[] contactContactivityForLoadStep26 = contactContactivityForEachStep[25].Values.ToArray();
            double[] contactContactivityForLoadStep27 = contactContactivityForEachStep[26].Values.ToArray();
            double[] contactContactivityForLoadStep28 = contactContactivityForEachStep[27].Values.ToArray();
            double[] contactContactivityForLoadStep29 = contactContactivityForEachStep[28].Values.ToArray();
            double[] contactContactivityForLoadStep30 = contactContactivityForEachStep[29].Values.ToArray();
            double[] contactContactivityForLoadStep31 = contactContactivityForEachStep[30].Values.ToArray();
            double[] contactContactivityForLoadStep32 = contactContactivityForEachStep[31].Values.ToArray();
            double[] contactContactivityForLoadStep33 = contactContactivityForEachStep[32].Values.ToArray();
            double[] contactContactivityForLoadStep34 = contactContactivityForEachStep[33].Values.ToArray();
            double[] contactContactivityForLoadStep35 = contactContactivityForEachStep[34].Values.ToArray();
            double[] contactContactivityForLoadStep36 = contactContactivityForEachStep[35].Values.ToArray();
            double[] contactContactivityForLoadStep37 = contactContactivityForEachStep[36].Values.ToArray();
            double[] contactContactivityForLoadStep38 = contactContactivityForEachStep[37].Values.ToArray();
            double[] contactContactivityForLoadStep39 = contactContactivityForEachStep[38].Values.ToArray();
            double[] contactContactivityForLoadStep40 = contactContactivityForEachStep[39].Values.ToArray();
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol1), fullThermalSolfinal1, @"C:\Users\Public\Documents\Results1.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol2), fullThermalSolfinal2, @"C:\Users\Public\Documents\Results2.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol3), fullThermalSolfinal3, @"C:\Users\Public\Documents\Results3.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol4), fullThermalSolfinal4, @"C:\Users\Public\Documents\Results4.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol5), fullThermalSolfinal5, @"C:\Users\Public\Documents\Results5.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol6), fullThermalSolfinal6, @"C:\Users\Public\Documents\Results6.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol7), fullThermalSolfinal7, @"C:\Users\Public\Documents\Results7.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol8), fullThermalSolfinal8, @"C:\Users\Public\Documents\Results8.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol9), fullThermalSolfinal9, @"C:\Users\Public\Documents\Results9.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol10), fullThermalSolfinal10, @"C:\Users\Public\Documents\Results10.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol11), fullThermalSolfinal11, @"C:\Users\Public\Documents\Results11.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol12), fullThermalSolfinal12, @"C:\Users\Public\Documents\Results12.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol13), fullThermalSolfinal13, @"C:\Users\Public\Documents\Results13.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol14), fullThermalSolfinal14, @"C:\Users\Public\Documents\Results14.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol15), fullThermalSolfinal15, @"C:\Users\Public\Documents\Results15.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol16), fullThermalSolfinal16, @"C:\Users\Public\Documents\Results16.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol17), fullThermalSolfinal17, @"C:\Users\Public\Documents\Results17.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol18), fullThermalSolfinal18, @"C:\Users\Public\Documents\Results18.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol19), fullThermalSolfinal19, @"C:\Users\Public\Documents\Results19.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol20), fullThermalSolfinal20, @"C:\Users\Public\Documents\Results20.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol21), fullThermalSolfinal21, @"C:\Users\Public\Documents\Results21.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol22), fullThermalSolfinal22, @"C:\Users\Public\Documents\Results22.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol23), fullThermalSolfinal23, @"C:\Users\Public\Documents\Results23.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol24), fullThermalSolfinal24, @"C:\Users\Public\Documents\Results24.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol25), fullThermalSolfinal25, @"C:\Users\Public\Documents\Results25.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol26), fullThermalSolfinal26, @"C:\Users\Public\Documents\Results26.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol27), fullThermalSolfinal27, @"C:\Users\Public\Documents\Results27.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol28), fullThermalSolfinal28, @"C:\Users\Public\Documents\Results28.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol29), fullThermalSolfinal29, @"C:\Users\Public\Documents\Results29.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol30), fullThermalSolfinal30, @"C:\Users\Public\Documents\Results30.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol31), fullThermalSolfinal31, @"C:\Users\Public\Documents\Results31.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol32), fullThermalSolfinal32, @"C:\Users\Public\Documents\Results32.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol33), fullThermalSolfinal33, @"C:\Users\Public\Documents\Results33.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol34), fullThermalSolfinal34, @"C:\Users\Public\Documents\Results34.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol35), fullThermalSolfinal35, @"C:\Users\Public\Documents\Results35.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol36), fullThermalSolfinal36, @"C:\Users\Public\Documents\Results36.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol37), fullThermalSolfinal37, @"C:\Users\Public\Documents\Results37.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol38), fullThermalSolfinal38, @"C:\Users\Public\Documents\Results38.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol39), fullThermalSolfinal39, @"C:\Users\Public\Documents\Results39.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol40), fullThermalSolfinal40, @"C:\Users\Public\Documents\Results40.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep1, @"C:\Users\Public\Documents\contactivity1.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep2, @"C:\Users\Public\Documents\contactivity2.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep3, @"C:\Users\Public\Documents\contactivity3.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep4, @"C:\Users\Public\Documents\contactivity4.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep5, @"C:\Users\Public\Documents\contactivity5.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep6, @"C:\Users\Public\Documents\contactivity6.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep7, @"C:\Users\Public\Documents\contactivity7.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep8, @"C:\Users\Public\Documents\contactivity8.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep9, @"C:\Users\Public\Documents\contactivity9.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep10, @"C:\Users\Public\Documents\contactivity10.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep11, @"C:\Users\Public\Documents\contactivity11.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep12, @"C:\Users\Public\Documents\contactivity12.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep13, @"C:\Users\Public\Documents\contactivity13.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep14, @"C:\Users\Public\Documents\contactivity14.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep15, @"C:\Users\Public\Documents\contactivity15.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep16, @"C:\Users\Public\Documents\contactivity16.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep17, @"C:\Users\Public\Documents\contactivity17.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep18, @"C:\Users\Public\Documents\contactivity18.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep19, @"C:\Users\Public\Documents\contactivity19.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep20, @"C:\Users\Public\Documents\contactivity20.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep21, @"C:\Users\Public\Documents\contactivity21.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep22, @"C:\Users\Public\Documents\contactivity22.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep23, @"C:\Users\Public\Documents\contactivity23.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep24, @"C:\Users\Public\Documents\contactivity24.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep25, @"C:\Users\Public\Documents\contactivity25.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep26, @"C:\Users\Public\Documents\contactivity26.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep27, @"C:\Users\Public\Documents\contactivity27.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep28, @"C:\Users\Public\Documents\contactivity28.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep29, @"C:\Users\Public\Documents\contactivity29.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep30, @"C:\Users\Public\Documents\contactivity30.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep31, @"C:\Users\Public\Documents\contactivity31.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep32, @"C:\Users\Public\Documents\contactivity32.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep33, @"C:\Users\Public\Documents\contactivity33.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep34, @"C:\Users\Public\Documents\contactivity34.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep35, @"C:\Users\Public\Documents\contactivity35.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep36, @"C:\Users\Public\Documents\contactivity36.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep37, @"C:\Users\Public\Documents\contactivity37.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep38, @"C:\Users\Public\Documents\contactivity38.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep39, @"C:\Users\Public\Documents\contactivity39.dat");
            VectorOperations.PrintVectorToFile(contactContactivityForLoadStep40, @"C:\Users\Public\Documents\contactivity40.dat");
            structuralSolutions.Add(fullStructuralSol1);
            structuralSolutions.Add(fullStructuralSol2);
            structuralSolutions.Add(fullStructuralSol3);
            structuralSolutions.Add(fullStructuralSol4);
            structuralSolutions.Add(fullStructuralSol5);
            structuralSolutions.Add(fullStructuralSol6);
            structuralSolutions.Add(fullStructuralSol7);
            structuralSolutions.Add(fullStructuralSol8);
            structuralSolutions.Add(fullStructuralSol9);
            structuralSolutions.Add(fullStructuralSol10);
            structuralSolutions.Add(fullStructuralSol11);
            structuralSolutions.Add(fullStructuralSol12);
            structuralSolutions.Add(fullStructuralSol13);
            structuralSolutions.Add(fullStructuralSol14);
            structuralSolutions.Add(fullStructuralSol15);
            structuralSolutions.Add(fullStructuralSol16);
            structuralSolutions.Add(fullStructuralSol17);
            structuralSolutions.Add(fullStructuralSol18);
            structuralSolutions.Add(fullStructuralSol19);
            structuralSolutions.Add(fullStructuralSol20);
            structuralSolutions.Add(fullStructuralSol21);
            structuralSolutions.Add(fullStructuralSol22);
            structuralSolutions.Add(fullStructuralSol23);
            structuralSolutions.Add(fullStructuralSol24);
            structuralSolutions.Add(fullStructuralSol25);
            structuralSolutions.Add(fullStructuralSol26);
            structuralSolutions.Add(fullStructuralSol27);
            structuralSolutions.Add(fullStructuralSol28);
            structuralSolutions.Add(fullStructuralSol29);
            structuralSolutions.Add(fullStructuralSol30);
            structuralSolutions.Add(fullStructuralSol31);
            structuralSolutions.Add(fullStructuralSol32);
            structuralSolutions.Add(fullStructuralSol33);
            structuralSolutions.Add(fullStructuralSol34);
            structuralSolutions.Add(fullStructuralSol35);
            structuralSolutions.Add(fullStructuralSol36);
            structuralSolutions.Add(fullStructuralSol37);
            structuralSolutions.Add(fullStructuralSol38);
            structuralSolutions.Add(fullStructuralSol39);
            structuralSolutions.Add(fullStructuralSol40);

            double[] Xvec1Final = new double[totalNodes / 2];
            double[] Yvec1Final = new double[totalNodes / 2];
            double[] Xvec2Final = new double[totalNodes / 2];
            double[] Yvec2Final = new double[totalNodes / 2];
            double[] Ζvec1Final = new double[totalNodes / 2];
            double[] Ζvec2Final = new double[totalNodes / 2];

            //Array.Copy(xFinalNodalCoor, 0, Xvec1Final, 0, totalNodes / 2);
            //Array.Copy(yFinalNodalCoor, 0, Yvec1Final, 0, totalNodes / 2);
            //Array.Copy(fullThermalSol4, 0, Ζvec1Final, 0, totalNodes / 2);
            //Array.Copy(xFinalNodalCoor, totalNodes / 2, Xvec2Final, 0, totalNodes / 2);
            //Array.Copy(yFinalNodalCoor, totalNodes / 2, Yvec2Final, 0, totalNodes / 2);
            //Array.Copy(fullThermalSol4, totalNodes / 2, Ζvec2Final, 0, totalNodes / 2);

            //List<HeatMapData> plots2 = new List<HeatMapData>();
            //plots2.Add(new HeatMapData() { Xcoordinates = Xvec1Final, Ycoordinates = Yvec1Final, Temperatures = Ζvec1Final });
            //plots2.Add(new HeatMapData() { Xcoordinates = Xvec2Final, Ycoordinates = Yvec2Final, Temperatures = Ζvec2Final });

            //ShowToGUI.PlotHeatMap(plots2);

            //string path = @"C:\Users\Public\Documents\Total\1final";
            //string path2 = @"C:\Users\Public\Documents\Total\2final";
            //ExportToFile.CreateContourDataForMatlab(Xvec1Final, Yvec1Final, Ζvec1Final, nodesInYCoor, nodesInXCoor, path);
            //ExportToFile.CreateContourDataForMatlab(Xvec2Final, Yvec2Final, Ζvec2Final, nodesInYCoor, nodesInXCoor, path2);

            //ExportToFile.ExportGeometryDataWithTemperatures(finalNodes, fullTempSol);

            GnuPlot.Close();

            //while (true)
            //{
            //    if (File.Exists(AppContext.BaseDirectory + "gnuplot.png") && new FileInfo(AppContext.BaseDirectory + "gnuplot.png").Length > 0)
            //    {
            //        break;

            //    }
            //    Thread.Sleep(100);
            //}
            //GnuPlot.KillProcess();
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

