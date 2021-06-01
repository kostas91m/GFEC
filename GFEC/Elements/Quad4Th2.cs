using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class Quad4Th2 : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        public double kc;
        private double A { get; set; }
        private double B { get; set; }
        //private double thickness = 1.0; //To be included in Element Properties
        //private double density = 1.0; //To be included in Element Properties

        public Quad4Th2(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, false, false, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, false, false, false, false, false };
            ElementFreedomSignature[3] = new bool[] { true, false, false, false, false, false };
            ElementFreedomSignature[4] = new bool[] { true, false, false, false, false, false };
            DisplacementVector = new double[4];
            A = properties.A;
            B = properties.B;
        }

        public double ClosestPointProjection()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }

        public Dictionary<int, INode> NodesAtFinalState()
        {
            throw new Exception("Method not implemenented");
        }
        public List<double[]> GetStressVector()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public List<double[]> GetStrainVector()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public List<double[]> GetGaussPointsInPhysicalSpace()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public List<double[]> GetStressFromElementsNodes()
        {
            throw new Exception("Method not implemenented");
        }
        public List<double[]> GetStrainFromElementsNodes()
        {
            throw new Exception("Method not implemenented");
        }
        public double[,] CreateGlobalStiffnessMatrix()
        {
            kc = Properties.ThermalConductivity;
            double[,] K = new double[4, 4];

            K[0, 0] = 2 * (Math.Pow(A, 2) + Math.Pow(B, 2)) * kc / (6.0 * A * B);
            K[0, 1] = (Math.Pow(A, 2) - 2 * Math.Pow(B, 2)) * kc / (6.0 * A * B);
            K[0, 2] = -(Math.Pow(A, 2) + Math.Pow(B, 2)) * kc / (6.0 * A * B);
            K[0, 3] = (Math.Pow(B, 2) - 2 * Math.Pow(A, 2)) * kc / (6.0 * A * B);

            K[1, 0] = K[0, 1];
            K[1, 1] = 2 * (Math.Pow(A, 2) + Math.Pow(B, 2)) * kc / (6.0 * A * B);
            K[1, 2] = (Math.Pow(B, 2) - 2 * Math.Pow(A, 2)) * kc / (6.0 * A * B);
            K[1, 3] = -(Math.Pow(A, 2) + Math.Pow(B, 2)) * kc / (6.0 * A * B);

            K[2, 0] = K[0, 2];
            K[2, 1] = K[1, 2];
            K[2, 2] = 2 * (Math.Pow(A, 2) + Math.Pow(B, 2)) * kc / (6.0 * A * B);
            K[2, 3] = (Math.Pow(A, 2) - 2 * Math.Pow(B, 2)) * kc / (6.0 * A * B);

            K[3, 0] = K[0, 3];
            K[3, 1] = K[1, 3];
            K[3, 2] = K[2, 3];
            K[3, 3] = 2 * (Math.Pow(A, 2) + Math.Pow(B, 2)) * kc / (6.0 * A * B);

            return K;
        }

        public double[,] CreateMassMatrix()
        {
            throw new Exception("Mass matrix not implemented for Quad4Th element");
        }

        public double[,] CreateDampingMatrix()
        {
            throw new Exception("Damping matrix not implemented for Quad4Th element");
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double[] intForces;
            double[,] stiff = CreateGlobalStiffnessMatrix();

            intForces = VectorOperations.MatrixVectorProduct(stiff, DisplacementVector);
            intForces = VectorOperations.VectorScalarProductNew(intForces, 1.0);
            return intForces;
        }
    }
}

