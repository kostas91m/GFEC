using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;

namespace GFEC
{
    public static class CoupledThermalStructuralCNTsInAngle6
    {
        private const int totalNodes = 486;
        private const int AddedNodes = 1;
        private const int RodElements = AddedNodes;
        private const int totalContactElements = 71;/*77;*//*71;*//*49;*//*30;*///10;//+11
        private const int totalElements = 320;
        private const int nodesInXCoor = 81;
        private const int nodesInYCoor = 3;
        private const double scaleFactor = 1.0;
        private const double xIntervals = 0.1;
        private const double yIntervals = 0.1;
        //offset for 10 contacts
        //private const double offset = 7.0 + 0.075;/*7.0;*//*7.0 + 0.05;*//*7.0 + 0.075;*//*7.0 + 0.10;*//*7.0 + 0.125;*//*7.0 + 0.15;*7.0 + 0.165*////- 0.05;//9 contacts
        private const double offset = 7.0 + 0.07 - 6.10; /*7.0 + 0.0695 - 6.10;*/ /*7.0 + 0.0695 - 15 * 0.0065 - 6.10;*//*7.0 + 0.0695;*//*7.0 + 0.07;*//*7.0 + 0.10;*//*7.0 + 0.11;*///- 6.10;//71 contacts
        //private const double offset = 7.0 + 0.070 - 6.10 - 0.60;/*7.0 + 0.0695 - 15 * 0.0065 - 6.10 - 0.60;*/ /*7.0 + 0.0695 - 6.10 - 0.60;*//*7.0 + 0.075 - 6.10 - 0.60;*//*7.0 + 0.09 - 6.10 - 0.60;*//*7.0 + 0.11 - 6.10 - 0.60;*//*7.0 + 0.135 - 6.10 - 0.60;*///-0.60//76 contacts
        //private const double offset = 7.0 + 0.070 - 3.90;/*7.0 + 0.0695 - 15 * 0.0065 - 3.90;*/ /*7.0 + 0.0695 - 3.90;*//*7.0 + 0.070 - 3.90;*//*7.0 + 0.075 - 3.90;*//*7.0 + 0.09 - 3.90;*//*7.0 + 0.11 - 3.90;*//*7.0 + 0.135 - 3.90;*///-3.90//48 contacts
        //private const double offset = 7.0 + 0.075 - 2;/*7.0 + 0.0695 - 15 * 0.0065 - 2;*/ /*7.0 + 0.0695 - 2;*//*7.0 + 0.070 - 2;*//*7.0 + 0.075 - 2;*//*7.0 + 0.09 - 2;*//*7.0 + 0.11 - 2;*//*7.0 + 0.135 - 2;*///-2//29 contacts
        //private const double offset = 7.0 + 0.0695 - 3.2;/*7.0 + 0.0695 - 15 * 0.0065 - 3.2;*/ /*7.0 + 0.0695 - 3.2;*//*7.0 + 0.070 - 3.2;*//*7.0 + 0.075 - 3.2;*//*7.0 + 0.09 - 3.2;*//*7.0 + 0.11 - 3.2;*//*7.0 + 0.135 - 3.2;*///-3.2//41 contacts
        //private const double gap = 1.14;
        public static ISolver structuralSolution;
        //private const double angle = Math.PI / 2.2;
        private const double angle = 0.491666667 * Math.PI - Math.PI * 0.016666667;/*Math.PI * 0.46666667 - 2 * Math.PI * 0.016666667;*//*Math.PI * 0.46666667 - Math.PI * 0.016666667;*//*0.491666667 * Math.PI - 2 * Math.PI * 0.016666667;*//*Math.PI * 0.46666667;*//* 0.491666667 * Math.PI - Math.PI * 0.016666667;*//*Math.PI * 0.48485;/* 0.49166667 * Math.PI;*/
        private const double gap = 0.628;/*1.6633;*//*1.252;*//*1.045;*//*0.837*//*0.628;*//*0.381*//*0.2095;*/
        private const int ThermalDof1 = 2;
        private const int ThermalDof2 = nodesInXCoor * (nodesInYCoor - 1) + 2;

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
        const double externalStructuralLoad = -200.0;//20 MPa
        const double T0 = 100.0;
        const double cond = 3300 * 1.0e-6;
        static double externalHeatLoad = -2 * T0 * (cond / (6 * xIntervals * yIntervals)) * ((Math.Pow(xIntervals, 2) - 2 * Math.Pow(yIntervals, 2)) - (Math.Pow(xIntervals, 2) + Math.Pow(yIntervals, 2)));
        //const double externalHeatLoad = 2500.0 * 1e-6;
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

            //for (int i = 0; i < nodesInYCoor; i++) //upper beam left side support
            //{
            //    boundedDofs.Add(i * nodesInXCoor * 2 + 1);
            //    boundedDofs.Add(i * nodesInXCoor * 2 + 2);
            //}

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

            for (int i = totalNodes + 1; i <= totalNodes + AddedNodes; i++)
            {
                boundedDofs.Add(i * 2);
                boundedDofs.Add(i * 2 - 1);
            }

            for (int i = totalNodes / 2 + 1; i <= totalNodes; i++)
            {
                boundedDofs.Add(i * 2 + 0);
                boundedDofs.Add(i * 2 - 1); //lower beam support for all nodes
            }
            boundedDofs.Add(82 * 2 - 1);


            structuralBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateThermalBoundaryConditions()
        {
            List<int> boundedDofs = new List<int>();
            for (int i = 0; i < nodesInYCoor; i++)
            {
                boundedDofs.Add(i * nodesInXCoor + 1); //upper beam left side support
            }
            for (int i = 0; i < nodesInYCoor; i++)
            {
                boundedDofs.Add(nodesInYCoor * nodesInXCoor + nodesInXCoor * (i + 1)); //lower beam right side support
            }
            thermalBoundaryConditions = boundedDofs.ToArray<int>();
        }

        private static void CreateStructuralLoadVector()
        {
            loadedStructuralDOFs = new List<int>();
            for (int i = 0; i < nodesInXCoor; i++)
            {
                loadedStructuralDOFs.Add(nodesInXCoor * nodesInYCoor * 2 - 2 * i);
            }
            externalForcesStructuralVector = new double[(totalNodes + AddedNodes) * 2];
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

            //extra nodes for rods on upper left side
            double ancos = Math.Cos(angle);
            nodes[487] = new Node(ancos * xIntervals * 2, 3 * yIntervals);
            //nodes[487] = new Node(-xIntervals, 0);
            //nodes[488] = new Node(-xIntervals, yIntervals);
            //nodes[489] = new Node(-xIntervals, 2 * yIntervals);
            //nodes[490] = new Node(0, -yIntervals);
            //nodes[491] = new Node(0.009515, 3 * yIntervals);
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

            //extra nodes for rods on upper left side
            //nodes[487] = new Node(0.009515, 3 * yIntervals);
            //nodes[487] = new Node(-xIntervals, 0);
            //nodes[488] = new Node(-xIntervals, yIntervals);
            //nodes[489] = new Node(-xIntervals, 2 * yIntervals);
            //nodes[490] = new Node(0, -yIntervals);
            //nodes[491] = new Node(0.009515, 3 * yIntervals);
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
            int localcounter = 0;
            for (int i = 1; i <= totalContactElements - 1; i++)
            {
                int lowerMiddleNode = new int();
                int lowerLeftNode = new int();
                int lowerRightNode = new int();
                int upperNode = new int();
                int lowerMiddleNode2 = new int();
                int lowerLeftNode2 = new int();
                int lowerRightNode2 = new int();
                int upperNode2 = new int();
                if (i == 1)
                {
                    lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i + 3;
                    lowerLeftNode = lowerMiddleNode - 1;
                    lowerRightNode = lowerMiddleNode + 1;
                    upperNode = nodesInXCoor - totalContactElements + i + 1;
                    localcounter += 1;
                    connectivity[totalElements + totalContactElements - 1 + localcounter] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
                }
                else if (i == totalContactElements - 1)
                {
                    lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i - 1;
                    lowerLeftNode = lowerMiddleNode - 1;
                    lowerRightNode = lowerMiddleNode + 1;
                    upperNode = nodesInXCoor - totalContactElements + i + 1;
                    localcounter += 1;
                    connectivity[totalElements + totalContactElements - 1 + localcounter] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
                }
                else if (i == 2)
                {
                    lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i;
                    lowerLeftNode = lowerMiddleNode - 1;
                    lowerRightNode = lowerMiddleNode;
                    upperNode = nodesInXCoor - totalContactElements + i + 1;
                    localcounter += 1;
                    connectivity[totalElements + totalContactElements - 1 + localcounter] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
                    lowerMiddleNode2 = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i + 3;
                    lowerLeftNode2 = lowerMiddleNode2 - 1;
                    lowerRightNode2 = lowerMiddleNode2 + 1;
                    upperNode2 = nodesInXCoor - totalContactElements + i + 1;
                    localcounter += 1;
                    connectivity[totalElements + totalContactElements - 1 + localcounter] = new Dictionary<int, int>() { { 1, lowerLeftNode2 }, { 2, lowerRightNode2 }, { 3, upperNode2 } };
                }
                else if(i == totalContactElements - 2)
                {
                    lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i - 1;
                    lowerLeftNode = lowerMiddleNode - 1;
                    lowerRightNode = lowerMiddleNode + 1;
                    upperNode = nodesInXCoor - totalContactElements + i + 1;
                    localcounter += 1;
                    connectivity[totalElements + totalContactElements - 1 + localcounter] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
                    lowerMiddleNode2 = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i + 2;
                    lowerLeftNode2 = lowerMiddleNode2 - 1;
                    lowerRightNode2 = lowerMiddleNode2;
                    upperNode2 = nodesInXCoor - totalContactElements + i + 1;
                    localcounter += 1;
                    connectivity[totalElements + totalContactElements - 1 + localcounter] = new Dictionary<int, int>() { { 1, lowerLeftNode2 }, { 2, lowerRightNode2 }, { 3, upperNode2 } };
                }
                else
                {
                    lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i - 1;
                    lowerLeftNode = lowerMiddleNode - 1;
                    lowerRightNode = lowerMiddleNode + 1;
                    upperNode = nodesInXCoor - totalContactElements + i + 1;
                    localcounter += 1;
                    connectivity[totalElements + totalContactElements - 1 + localcounter] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
                    lowerMiddleNode2 = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i + 3;
                    lowerLeftNode2 = lowerMiddleNode2 - 1;
                    lowerRightNode2 = lowerMiddleNode2 + 1;
                    upperNode2 = nodesInXCoor - totalContactElements + i + 1;
                    localcounter += 1;
                    connectivity[totalElements + totalContactElements - 1 + localcounter] = new Dictionary<int, int>() { { 1, lowerLeftNode2 }, { 2, lowerRightNode2 }, { 3, upperNode2 } };
                }

            }

            //Rod elements
            int count = connectivity.Count;
            connectivity[count + 1] = new Dictionary<int, int>() { { 1, 163 }, { 2, totalNodes + AddedNodes } };
            //connectivity[count + 1] = new Dictionary<int, int>() { { 1, 1 }, { 2, 487 } };
            //connectivity[count + 2] = new Dictionary<int, int>() { { 1, 82 }, { 2, 488 } };
            //connectivity[count + 3] = new Dictionary<int, int>() { { 1, 163 }, { 2, 489 } };
            //connectivity[count + 4] = new Dictionary<int, int>() { { 1, 1 }, { 2, 490 } };
            //connectivity[count + 5] = new Dictionary<int, int>() { { 1, 163 }, { 2, 491 } };

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
            for (int i = 1; i <= totalContactElements - 1; i++)
            {
                int lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i + 1;
                int lowerLeftNode = lowerMiddleNode - 1;
                int lowerRightNode = lowerMiddleNode + 1;
                int upperNode = nodesInXCoor - totalContactElements + i + 1;
                connectivity[totalElements + i] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
            }
            int localcounter = 0;
            for (int i = 1; i <= totalContactElements - 1; i++)
            {
                int lowerMiddleNode = new int();
                int lowerLeftNode = new int();
                int lowerRightNode = new int();
                int upperNode = new int();
                int lowerMiddleNode2 = new int();
                int lowerLeftNode2 = new int();
                int lowerRightNode2 = new int();
                int upperNode2 = new int();
                if (i == 1)
                {
                    lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i + 3;
                    lowerLeftNode = lowerMiddleNode - 1;
                    lowerRightNode = lowerMiddleNode + 1;
                    upperNode = nodesInXCoor - totalContactElements + i + 1;
                    localcounter += 1;
                    connectivity[totalElements + totalContactElements - 1 + localcounter] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
                }
                else if (i == totalContactElements - 1)
                {
                    lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i - 1;
                    lowerLeftNode = lowerMiddleNode - 1;
                    lowerRightNode = lowerMiddleNode + 1;
                    upperNode = nodesInXCoor - totalContactElements + i + 1;
                    localcounter += 1;
                    connectivity[totalElements + totalContactElements - 1 + localcounter] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
                }
                else
                {
                    lowerMiddleNode = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i - 1;
                    lowerLeftNode = lowerMiddleNode - 1;
                    lowerRightNode = lowerMiddleNode + 1;
                    upperNode = nodesInXCoor - totalContactElements + i + 1;
                    localcounter += 1;
                    connectivity[totalElements + totalContactElements - 1 + localcounter] = new Dictionary<int, int>() { { 1, lowerLeftNode }, { 2, lowerRightNode }, { 3, upperNode } };
                    lowerMiddleNode2 = 2 * nodesInXCoor * nodesInYCoor - nodesInXCoor + i + 3;
                    lowerLeftNode2 = lowerMiddleNode - 1;
                    lowerRightNode2 = lowerMiddleNode + 1;
                    upperNode2 = nodesInXCoor - totalContactElements + i + 1;
                    localcounter += 1;
                    connectivity[totalElements + totalContactElements - 1 + localcounter] = new Dictionary<int, int>() { { 1, lowerLeftNode2 }, { 2, lowerRightNode2 }, { 3, upperNode2 } };
                }
            }
            //Rod elements
            //int count = connectivity.Count;
            //connectivity[count + 1] = new Dictionary<int, int>() { { 1, 163 }, { 2, totalNodes + AddedNodes } };
            //connectivity[count + 1] = new Dictionary<int, int>() { { 1, 1 }, { 2, 487 } };
            //connectivity[count + 2] = new Dictionary<int, int>() { { 1, 82 }, { 2, 488 } };
            //connectivity[count + 3] = new Dictionary<int, int>() { { 1, 163 }, { 2, 489 } };
            //connectivity[count + 4] = new Dictionary<int, int>() { { 1, 1 }, { 2, 490 } };
            //connectivity[count + 5] = new Dictionary<int, int>() { { 1, 163 }, { 2, 491 } };

            return connectivity;
        }

        private static Dictionary<int, bool[]> CreateNodeFAT()
        {
            Dictionary<int, bool[]> nodeFAT = new Dictionary<int, bool[]>();
            for (int i = 1; i <= totalNodes + AddedNodes; i++)
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

            for (int i = totalElements + 1; i <= totalElements + 3 * totalContactElements - 5; i++)
            {
                elementProperties[i] = new ElementProperties(E, A, type2);
                elementProperties[i].Density = density;
                elementProperties[i].Thickness = thickness;
            }

            int count = elementProperties.Count;
            for (int i = count + 1; i <= count + AddedNodes; i++)
            {
                elementProperties[i] = new ElementProperties(E / 5000, A, type3);
            }
            //elementProperties[count + 4] = new ElementProperties(E / 10000, A, type3);
            //elementProperties[count + 5] = new ElementProperties(E / 10000, A, type3);

            return elementProperties;
        }

        private static Dictionary<int, IElementProperties> CreateThermalElementProperties()
        {
            double thermalCond = solidThermalCond;
            double A = area;
            string type = "Quad4Th";
            string type2 = "ContactNtS2DTh";
            //string type3 = "ContactNtN2DTh";

            Dictionary<int, IElementProperties> elementProperties = new Dictionary<int, IElementProperties>();
            for (int i = 1; i <= totalElements; i++)
            {
                elementProperties[i] = new ElementProperties();
                elementProperties[i].ElementType = type;
                elementProperties[i].ThermalConductivity = thermalCond;
            }
            for (int i = totalElements + 1; i <= totalElements + 3 * totalContactElements - 5; i++)
            {
                elementProperties[i] = new ElementProperties();
                elementProperties[i].ElementType = type2;
                elementProperties[i].ThermalConductivity = thermalCond;
                elementProperties[i].SectionArea = A;
                elementProperties[i].SurfaceRoughness = roughness;
                elementProperties[i].ContactThermalConductivity = contactCond;
                elementProperties[i].YieldStrength = yieldStrength;
            }
            //for (int i = totalElements + totalContactElements; i < totalElements + totalContactElements + RodElements; i++)
            //{
            //    elementProperties[i] = new ElementProperties();
            //    elementProperties[i].ElementType = type3;
            //    elementProperties[i].ThermalConductivity = 0;
            //    elementProperties[i].SectionArea = A;
            //    elementProperties[i].SurfaceRoughness = roughness;
            //    elementProperties[i].ContactThermalConductivity = contactCond;
            //    elementProperties[i].YieldStrength = yieldStrength;
            //}
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
            structuralSolution.LinearScheme = new PCGSolver();
            //structuralSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
            structuralSolution.NonLinearScheme.Tolerance = 1e-5;
            structuralSolution.ActivateNonLinearSolver = true;
            structuralSolution.NonLinearScheme.numberOfLoadSteps = 120;

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
                for (int j = totalElements + 1; j <= totalElements + 3 * totalContactElements - 5; j++)
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
            List<double[]> thermalSolutions = new List<double[]>();
            //Dictionary<int, double[]> thermalSolutions = new Dictionary<int, double[]>();
            List<Dictionary<int, double>> contactContactivityForEachStep = new List<Dictionary<int, double>>();
            for (int k = 1; k <= allStepsSolutions.Count; k++)
            {
                IAssembly elementsAssembly2 = CreateThermalAssembly();

                for (int j = totalElements + 1; j <= totalElements + 3 * totalContactElements - 5; j++)
                {
                    double[] contactForce = allStepsContactForces[k][j];
                    elementsAssembly2.ElementsProperties[j].ContactForceValue = -contactForce[5];
                    double projectionPoint = allStepsProjectionPoints[k][j];
                    elementsAssembly2.ElementsProperties[j].Dx1 = projectionPoint;
                }

                elementsAssembly2.CreateElementsAssembly();
                elementsAssembly2.ActivateBoundaryConditions = true;
                double[,] globalStiffnessMatrix2 = elementsAssembly2.CreateTotalStiffnessMatrix();

                ISolver thermalSolution = new StaticSolver();
                thermalSolution.LinearScheme = new LUFactorization();
                thermalSolution.NonLinearScheme = new LoadControlledNewtonRaphson();
                thermalSolution.NonLinearScheme.Tolerance = 1e-9;
                thermalSolution.ActivateNonLinearSolver = true;
                thermalSolution.NonLinearScheme.numberOfLoadSteps = 40;

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
                //foreach (var dof in loadedThermalDOFs)
                //{
                //    externalHeatFlux[dof - 1] = externalHeatLoad;
                //}
                //for (int i = 61; i <= 75; i++)
                //{
                //    externalHeatFlux[61] = externalHeatLoad;
                //}
                double[] reducedExternalHeatFlux = BoundaryConditionsImposition.ReducedVector(externalHeatFlux, thermalSolution.AssemblyData.BoundedDOFsVector);
                thermalSolution.Solve(reducedExternalHeatFlux);
                double[] tempSol = thermalSolution.GetSolution();
                thermalSolutions.Add(tempSol);
                //thermalSolutions.Add(k, tempSol);

                Dictionary<int, double> contactContactivity = AssemblyHelpMethods.RetrieveContactContactivity(thermalSolution.AssemblyData);
                contactContactivityForEachStep.Add(contactContactivity);
            }
            //ExportToFile.ExportGeometryDataWithTemperatures(structuralSolution, thermalSolutions, thermalBoundaryConditions, @"C:\Users\Public\Documents\");
            //ExportToFile.ExportCondactivityForAllLoadSteps(contactContactivityForEachStep);

            int[] thermalBoundCond = thermalBoundaryConditions;
            double[] fullStructuralSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[3], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[6], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[9], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[12], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[15], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol6 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[18], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol7 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[21], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol8 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[24], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol9 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[27], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol10 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[30], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol11 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[33], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol12 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[36], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol13 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[39], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol14 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[42], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol15 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[45], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol16 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[48], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol17 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[51], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol18 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[54], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol19 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[57], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol20 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[60], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol21 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[63], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol22 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[66], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol23 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[69], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol24 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[72], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol25 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[75], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol26 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[78], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol27 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[81], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol28 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[84], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol29 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[87], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol30 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[90], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol31 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[93], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol32 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[96], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol33 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[99], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol34 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[102], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol35 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[105], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol36 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[108], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol37 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[111], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol38 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[114], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol39 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[117], elementsAssembly.BoundedDOFsVector);
            double[] fullStructuralSol40 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(allStepsSolutions[120], elementsAssembly.BoundedDOFsVector);
            double[] fullThermalSol1 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[2], thermalBoundCond);
            double[] fullThermalSol2 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[5], thermalBoundCond);
            double[] fullThermalSol3 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[8], thermalBoundCond);
            double[] fullThermalSol4 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[11], thermalBoundCond);
            double[] fullThermalSol5 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[14], thermalBoundCond);
            double[] fullThermalSol6 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[17], thermalBoundCond);
            double[] fullThermalSol7 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[20], thermalBoundCond);
            double[] fullThermalSol8 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[23], thermalBoundCond);
            double[] fullThermalSol9 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[26], thermalBoundCond);
            double[] fullThermalSol10 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[29], thermalBoundCond);
            double[] fullThermalSol11 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[32], thermalBoundCond);
            double[] fullThermalSol12 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[35], thermalBoundCond);
            double[] fullThermalSol13 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[38], thermalBoundCond);
            double[] fullThermalSol14 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[41], thermalBoundCond);
            double[] fullThermalSol15 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[44], thermalBoundCond);
            double[] fullThermalSol16 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[47], thermalBoundCond);
            double[] fullThermalSol17 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[50], thermalBoundCond);
            double[] fullThermalSol18 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[53], thermalBoundCond);
            double[] fullThermalSol19 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[56], thermalBoundCond);
            double[] fullThermalSol20 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[59], thermalBoundCond);
            double[] fullThermalSol21 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[62], thermalBoundCond);
            double[] fullThermalSol22 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[65], thermalBoundCond);
            double[] fullThermalSol23 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[68], thermalBoundCond);
            double[] fullThermalSol24 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[71], thermalBoundCond);
            double[] fullThermalSol25 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[74], thermalBoundCond);
            double[] fullThermalSol26 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[77], thermalBoundCond);
            double[] fullThermalSol27 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[80], thermalBoundCond);
            double[] fullThermalSol28 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[83], thermalBoundCond);
            double[] fullThermalSol29 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[86], thermalBoundCond);
            double[] fullThermalSol30 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[89], thermalBoundCond);
            double[] fullThermalSol31 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[92], thermalBoundCond);
            double[] fullThermalSol32 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[95], thermalBoundCond);
            double[] fullThermalSol33 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[98], thermalBoundCond);
            double[] fullThermalSol34 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[101], thermalBoundCond);
            double[] fullThermalSol35 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[104], thermalBoundCond);
            double[] fullThermalSol36 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[107], thermalBoundCond);
            double[] fullThermalSol37 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[110], thermalBoundCond);
            double[] fullThermalSol38 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[113], thermalBoundCond);
            double[] fullThermalSol39 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[116], thermalBoundCond);
            double[] fullThermalSol40 = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(thermalSolutions[119], thermalBoundCond);
            double[] fullThermalSolfinal1 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal2 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal3 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal4 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal5 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal6 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal7 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal8 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal9 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal10 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal11 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal12 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal13 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal14 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal15 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal16 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal17 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal18 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal19 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal20 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal21 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal22 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal23 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal24 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal25 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal26 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal27 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal28 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal29 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal30 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal31 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal32 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal33 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal34 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal35 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal36 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal37 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal38 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal39 = new double[fullThermalSol1.Length + 1];
            double[] fullThermalSolfinal40 = new double[fullThermalSol1.Length + 1];

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
            double[] contactContactivityForLoadStep1 = contactContactivityForEachStep[1].Values.ToArray();
            double[] contactContactivityForLoadStep2 = contactContactivityForEachStep[3].Values.ToArray();
            double[] contactContactivityForLoadStep3 = contactContactivityForEachStep[5].Values.ToArray();
            double[] contactContactivityForLoadStep4 = contactContactivityForEachStep[7].Values.ToArray();
            double[] contactContactivityForLoadStep5 = contactContactivityForEachStep[9].Values.ToArray();
            double[] contactContactivityForLoadStep6 = contactContactivityForEachStep[11].Values.ToArray();
            double[] contactContactivityForLoadStep7 = contactContactivityForEachStep[13].Values.ToArray();
            double[] contactContactivityForLoadStep8 = contactContactivityForEachStep[15].Values.ToArray();
            double[] contactContactivityForLoadStep9 = contactContactivityForEachStep[17].Values.ToArray();
            double[] contactContactivityForLoadStep10 = contactContactivityForEachStep[19].Values.ToArray();
            double[] contactContactivityForLoadStep11 = contactContactivityForEachStep[21].Values.ToArray();
            double[] contactContactivityForLoadStep12 = contactContactivityForEachStep[23].Values.ToArray();
            double[] contactContactivityForLoadStep13 = contactContactivityForEachStep[25].Values.ToArray();
            double[] contactContactivityForLoadStep14 = contactContactivityForEachStep[27].Values.ToArray();
            double[] contactContactivityForLoadStep15 = contactContactivityForEachStep[29].Values.ToArray();
            double[] contactContactivityForLoadStep16 = contactContactivityForEachStep[31].Values.ToArray();
            double[] contactContactivityForLoadStep17 = contactContactivityForEachStep[33].Values.ToArray();
            double[] contactContactivityForLoadStep18 = contactContactivityForEachStep[35].Values.ToArray();
            double[] contactContactivityForLoadStep19 = contactContactivityForEachStep[37].Values.ToArray();
            double[] contactContactivityForLoadStep20 = contactContactivityForEachStep[39].Values.ToArray();
            double[] contactContactivityForLoadStep21 = contactContactivityForEachStep[41].Values.ToArray();
            double[] contactContactivityForLoadStep22 = contactContactivityForEachStep[43].Values.ToArray();
            double[] contactContactivityForLoadStep23 = contactContactivityForEachStep[45].Values.ToArray();
            double[] contactContactivityForLoadStep24 = contactContactivityForEachStep[47].Values.ToArray();
            double[] contactContactivityForLoadStep25 = contactContactivityForEachStep[49].Values.ToArray();
            double[] contactContactivityForLoadStep26 = contactContactivityForEachStep[51].Values.ToArray();
            double[] contactContactivityForLoadStep27 = contactContactivityForEachStep[53].Values.ToArray();
            double[] contactContactivityForLoadStep28 = contactContactivityForEachStep[55].Values.ToArray();
            double[] contactContactivityForLoadStep29 = contactContactivityForEachStep[57].Values.ToArray();
            double[] contactContactivityForLoadStep30 = contactContactivityForEachStep[59].Values.ToArray();
            double[] contactContactivityForLoadStep31 = contactContactivityForEachStep[61].Values.ToArray();
            double[] contactContactivityForLoadStep32 = contactContactivityForEachStep[63].Values.ToArray();
            double[] contactContactivityForLoadStep33 = contactContactivityForEachStep[65].Values.ToArray();
            double[] contactContactivityForLoadStep34 = contactContactivityForEachStep[67].Values.ToArray();
            double[] contactContactivityForLoadStep35 = contactContactivityForEachStep[69].Values.ToArray();
            double[] contactContactivityForLoadStep36 = contactContactivityForEachStep[71].Values.ToArray();
            double[] contactContactivityForLoadStep37 = contactContactivityForEachStep[73].Values.ToArray();
            double[] contactContactivityForLoadStep38 = contactContactivityForEachStep[75].Values.ToArray();
            double[] contactContactivityForLoadStep39 = contactContactivityForEachStep[77].Values.ToArray();
            double[] contactContactivityForLoadStep40 = contactContactivityForEachStep[79].Values.ToArray();
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol1), fullThermalSolfinal1, @"C:\Users\Public\Documents\Yovanovitch1.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol2), fullThermalSolfinal2, @"C:\Users\Public\Documents\Yovanovitch2.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol3), fullThermalSolfinal3, @"C:\Users\Public\Documents\Yovanovitch3.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol4), fullThermalSolfinal4, @"C:\Users\Public\Documents\Yovanovitch4.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol5), fullThermalSolfinal5, @"C:\Users\Public\Documents\Yovanovitch5.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol6), fullThermalSolfinal6, @"C:\Users\Public\Documents\Yovanovitch6.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol7), fullThermalSolfinal7, @"C:\Users\Public\Documents\Yovanovitch7.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol8), fullThermalSolfinal8, @"C:\Users\Public\Documents\Yovanovitch8.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol9), fullThermalSolfinal9, @"C:\Users\Public\Documents\Yovanovitch9.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol10), fullThermalSolfinal10, @"C:\Users\Public\Documents\Yovanovitch10.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol11), fullThermalSolfinal11, @"C:\Users\Public\Documents\Yovanovitch11.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol12), fullThermalSolfinal12, @"C:\Users\Public\Documents\Yovanovitch12.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol13), fullThermalSolfinal13, @"C:\Users\Public\Documents\Yovanovitch13.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol14), fullThermalSolfinal14, @"C:\Users\Public\Documents\Yovanovitch14.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol15), fullThermalSolfinal15, @"C:\Users\Public\Documents\Yovanovitch15.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol16), fullThermalSolfinal16, @"C:\Users\Public\Documents\Yovanovitch16.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol17), fullThermalSolfinal17, @"C:\Users\Public\Documents\Yovanovitch17.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol18), fullThermalSolfinal18, @"C:\Users\Public\Documents\Yovanovitch18.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol19), fullThermalSolfinal19, @"C:\Users\Public\Documents\Yovanovitch19.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol20), fullThermalSolfinal20, @"C:\Users\Public\Documents\Yovanovitch20.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol21), fullThermalSolfinal21, @"C:\Users\Public\Documents\Yovanovitch21.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol22), fullThermalSolfinal22, @"C:\Users\Public\Documents\Yovanovitch22.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol23), fullThermalSolfinal23, @"C:\Users\Public\Documents\Yovanovitch23.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol24), fullThermalSolfinal24, @"C:\Users\Public\Documents\Yovanovitch24.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol25), fullThermalSolfinal25, @"C:\Users\Public\Documents\Yovanovitch25.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol26), fullThermalSolfinal26, @"C:\Users\Public\Documents\Yovanovitch26.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol27), fullThermalSolfinal27, @"C:\Users\Public\Documents\Yovanovitch27.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol28), fullThermalSolfinal28, @"C:\Users\Public\Documents\Yovanovitch28.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol29), fullThermalSolfinal29, @"C:\Users\Public\Documents\Yovanovitch29.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol30), fullThermalSolfinal30, @"C:\Users\Public\Documents\Yovanovitch30.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol31), fullThermalSolfinal31, @"C:\Users\Public\Documents\Yovanovitch31.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol32), fullThermalSolfinal32, @"C:\Users\Public\Documents\Yovanovitch32.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol33), fullThermalSolfinal33, @"C:\Users\Public\Documents\Yovanovitch33.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol34), fullThermalSolfinal34, @"C:\Users\Public\Documents\Yovanovitch34.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol35), fullThermalSolfinal35, @"C:\Users\Public\Documents\Yovanovitch35.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol36), fullThermalSolfinal36, @"C:\Users\Public\Documents\Yovanovitch36.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol37), fullThermalSolfinal37, @"C:\Users\Public\Documents\Yovanovitch37.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol38), fullThermalSolfinal38, @"C:\Users\Public\Documents\Yovanovitch38.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol39), fullThermalSolfinal39, @"C:\Users\Public\Documents\Yovanovitch39.dat");
            ExportToFile.ExportGeometryDataWithTemperatures(Assembly.CalculateFinalNodalCoordinates(elementsAssembly.Nodes, fullStructuralSol40), fullThermalSolfinal40, @"C:\Users\Public\Documents\Yovanovitch40.dat");
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

