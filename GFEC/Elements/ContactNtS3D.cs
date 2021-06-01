using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    class ContactNtS3D : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        private double PenaltyFactor { get; set; }
        double[] lastKsiVector;

        public ContactNtS3D(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[3] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[4] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[5] = new bool[] { true, true, true, false, false, false };
            DisplacementVector = new double[15];
            PenaltyFactor = properties.YoungMod * 1.0;
            lastKsiVector = new double[2];
        }

        public Dictionary<int, INode> NodesAtFinalState()
        {
            Dictionary<int, INode> finalNodes = new Dictionary<int, INode>();
            finalNodes[1] = new Node(Nodes[1].XCoordinate + DisplacementVector[0], Nodes[1].YCoordinate + DisplacementVector[1], Nodes[1].ZCoordinate + DisplacementVector[2]);
            finalNodes[2] = new Node(Nodes[2].XCoordinate + DisplacementVector[3], Nodes[2].YCoordinate + DisplacementVector[4], Nodes[2].ZCoordinate + DisplacementVector[5]);
            finalNodes[3] = new Node(Nodes[3].XCoordinate + DisplacementVector[6], Nodes[3].YCoordinate + DisplacementVector[7], Nodes[3].ZCoordinate + DisplacementVector[8]);
            finalNodes[4] = new Node(Nodes[4].XCoordinate + DisplacementVector[9], Nodes[4].YCoordinate + DisplacementVector[10], Nodes[4].ZCoordinate + DisplacementVector[11]);
            finalNodes[5] = new Node(Nodes[5].XCoordinate + DisplacementVector[12], Nodes[5].YCoordinate + DisplacementVector[13], Nodes[5].ZCoordinate + DisplacementVector[14]);
            return finalNodes;
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
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0, 0.0 });
            return l;
        }
        public List<double[]> GetStressFromElementsNodes()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0, 0.0 });
            return l;
        }
        public List<double[]> GetStrainFromElementsNodes()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0, 0.0 });
            return l;
        }
        private double[] xUpdatedVector()
        {
            Dictionary<int, INode> finalNodes = NodesAtFinalState();
            double[] xVectorUpdated = new double[15];
            for (int i = 0; i < 5; i++)
            {
                xVectorUpdated[3 * i] = finalNodes[i + 1].XCoordinate;
                xVectorUpdated[3 * i + 1] = finalNodes[i + 1].YCoordinate;
                xVectorUpdated[3 * i + 2] = finalNodes[i + 1].ZCoordinate;
            }
            return xVectorUpdated;
        }

        private Tuple<double[,], double[,], double[,]> CalculatePositionMatrix(double ksi1, double ksi2)
        {
            double N1 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 - ksi2);
            double N2 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 - ksi2);
            double N3 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 + ksi2);
            double N4 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 + ksi2);

            double dN11 = -1.0 / 4.0 * (1.0 - ksi2);
            double dN21 = 1.0 / 4.0 * (1.0 - ksi2);
            double dN31 = 1.0 / 4.0 * (1.0 + ksi2);
            double dN41 = -1.0 / 4.0 * (1.0 + ksi2);

            double dN12 = -1.0 / 4.0 * (1.0 - ksi1);
            double dN22 = -1.0 / 4.0 * (1.0 + ksi1);
            double dN32 = 1.0 / 4.0 * (1.0 + ksi1);
            double dN42 = 1.0 / 4.0 * (1.0 - ksi1);

            double[,] aMatrix = new double[,]
                {
                    { -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, 1.0, 0.0, 0.0 },
                    { 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, 1.0, 0.0 },
                    { 0.0, 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, 1.0 }
                };

            double[,] da1Matrix = new double[,]
                {
                    { -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 }
                };

            double[,] da2Matrix = new double[,]
                {
                    { -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 }
                };
            return new Tuple<double[,], double[,], double[,]>(aMatrix, da1Matrix, da2Matrix);
        }

        private List<double[]> SurfaceVectors(double ksi1, double ksi2)
        {
            Tuple<double[,], double[,], double[,]> positionsMatrices = CalculatePositionMatrix(ksi1, ksi2);
            double[,] da1 = positionsMatrices.Item2;
            double[,] da2 = positionsMatrices.Item3;
            double[] xUpdated = xUpdatedVector();

            List<double[]> dRho = new List<double[]>();
            dRho.Add(VectorOperations.VectorScalarProductNew(
                VectorOperations.MatrixVectorProduct(da1, xUpdated), -1.0));
            dRho.Add(VectorOperations.VectorScalarProductNew(
                VectorOperations.MatrixVectorProduct(da2, xUpdated), -1.0));
            return dRho;
        }

        private double[,] MetricTensor(List<double[]> dRho)
        {
            double[,] m = new double[2, 2];
            m[0, 0] = VectorOperations.VectorDotProduct(dRho[0], dRho[0]);
            m[0, 1] = VectorOperations.VectorDotProduct(dRho[0], dRho[1]);
            m[1, 0] = VectorOperations.VectorDotProduct(dRho[1], dRho[0]);
            m[1, 1] = VectorOperations.VectorDotProduct(dRho[1], dRho[1]);
            return m;
        }

        private double MetricTensorDet(double[,] m)
        {
            double detm = m[0, 0] * m[1, 1] - m[1, 0] * m[0, 1];

            return detm;
        }

        private double[,] InverseMetricTensor(double[,] m)
        {
            double detm = MetricTensorDet(m);

            double[,] mInv = MatrixOperations.ScalarMatrixProductNew(1.0 / detm,
                new double[,] {
                    { m[1, 1], -m[0, 1] },
                    {-m[1,0], m[0,0] }
                });
            return mInv;
        }
        private double[] NormalVector(double[,] m, List<double[]> dRho)
        {
            double detm = MetricTensorDet(m);
            double[] drhoBydrho = VectorOperations.VectorCrossProduct(dRho[0], dRho[1]);
            double[] n = VectorOperations.VectorScalarProductNew(drhoBydrho, 1.0 / (Math.Sqrt(detm)));
            return n;
        }

        private double CalculatePenetration(double ksi1, double ksi2)
        {
            double[] xUpdated = xUpdatedVector();
            List<double[]> dRho = SurfaceVectors(ksi1, ksi2);
            double[,] m = MetricTensor(dRho);
            double[] n = NormalVector(m, dRho);
            Tuple<double[,], double[,], double[,]> aMatrices = CalculatePositionMatrix(ksi1, ksi2);
            double[,] a = aMatrices.Item1;
            double[,] aT = MatrixOperations.Transpose(a);
            double ksi3 = VectorOperations.VectorDotProduct(
                xUpdated, VectorOperations.MatrixVectorProduct(
                    aT, n));
            return ksi3;
        }

        private double[] Calculate_f(List<double[]> dRho, double[,] aMatrix, double[] xUpdated)
        {
            double[] f = new double[2];
            f[0] = VectorOperations.VectorDotProduct(
                dRho[0], VectorOperations.MatrixVectorProduct(
                    aMatrix, xUpdated));
            f[1] = VectorOperations.VectorDotProduct(
                dRho[1], VectorOperations.MatrixVectorProduct(
                    aMatrix, xUpdated));
            return f;
        }

        private double Calculate_e(double[,] aMatrix, double[] xUpdated)
        {

            double[] dRho12 = new double[] {
                0.25*(xUpdated[0] - xUpdated[3] + xUpdated[6] - xUpdated[9]),
                0.25*(xUpdated[1] - xUpdated[4] + xUpdated[7] - xUpdated[10]),
                0.25*(xUpdated[2] - xUpdated[5] + xUpdated[8] - xUpdated[11])
            };

            double e = VectorOperations.VectorDotProduct(
                dRho12, VectorOperations.MatrixVectorProduct(
                    aMatrix, xUpdated));
            return e;
        }

        private double[] CalculateDeltaKsi(double detm, double[,] mTensor, double[] fVector, double e)
        {
            double scalar = 1.0 / (detm - Math.Pow(e, 2) + 2.0 * e * mTensor[0, 1]);
            double[,] matrix = new double[,]
            {
                {mTensor[1,1], e-mTensor[0,1] },
                {e-mTensor[1,0], mTensor[0,0] }
            };

            double[] deltaKsi = VectorOperations.VectorScalarProductNew(
                VectorOperations.MatrixVectorProduct(matrix, fVector),
                scalar);

            return deltaKsi;
        }

        private double[] Project(double[] ksiVectorInitial)
        {
            int maxIterations = 1000;
            double tol = Math.Pow(10.0, -4.0);
            double[] deltaKsi = new double[2];
            double norm = new double();
            double[] ksiVector = ksiVectorInitial;
            double[] xUpdated = xUpdatedVector();
            for(int i = 1; i <= maxIterations; i++)
            {
                double[] oldksiVector = ksiVector;
                Tuple<double[,], double[,], double[,]> aMatrices = CalculatePositionMatrix(ksiVector[0], ksiVector[1]);
                List<double[]> dRho = SurfaceVectors(ksiVector[0], ksiVector[1]);
                double[] f = Calculate_f(dRho, aMatrices.Item1, xUpdated);
                double e = Calculate_e(aMatrices.Item1, xUpdated);//double e = Calculate_e(aMatrices.Item2, xUpdated);
                double[,] m = MetricTensor(dRho);
                double detm = MetricTensorDet(m);
                deltaKsi = CalculateDeltaKsi(detm, m, f, e);
                ksiVector = VectorOperations.VectorVectorAddition(ksiVector, deltaKsi);//Iterations for the CPP?
                norm = VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(ksiVector, oldksiVector));
                if(norm <= tol)
                {
                    break;
                }
            }
            if (norm > tol)
            {
                throw new Exception("CPP not found in current iterations");
            }
            else
            {
                return ksiVector;

            }
        }

        private double CalculatePenetration(double[] normalVector, double[,] aMatrix, double[] xUpdated)
        {
            double ksi3 = VectorOperations.VectorDotProduct(
                xUpdated, VectorOperations.MatrixVectorProduct(
                    MatrixOperations.Transpose(aMatrix), normalVector));
            return ksi3;
        }

        private double[,] MainStiffnessPart(double penaltyFactor, double[] normalVector, double[,] aMatrix)
        {
            double[,] nxn = VectorOperations.VectorVectorTensorProduct(normalVector, normalVector);
            double[,] aT = MatrixOperations.Transpose(aMatrix);
            double[,] nxna = MatrixOperations.MatrixProduct(nxn, aMatrix);
            double[,] aTnxna = MatrixOperations.MatrixProduct(aT, nxna);
            double[,] Kmain = MatrixOperations.ScalarMatrixProductNew(penaltyFactor, aTnxna);
            return Kmain;
        }
        private double[,] RotationalStiffnessPart(double penaltyFactor, double[] normalVector, double[,] aMatrix, double[,] a1Matrix, double[,] a2Matrix, List<double[]> dRho, double ksi3)
        {
            double[,] m = MetricTensor(dRho);
            double[,] mInv = InverseMetricTensor(m);

            double scalar1 = penaltyFactor * ksi3 * mInv[0, 0];
            double scalar2 = penaltyFactor * ksi3 * mInv[1, 0];
            double scalar3 = penaltyFactor * ksi3 * mInv[0, 1];
            double scalar4 = penaltyFactor * ksi3 * mInv[1, 1];

            double[,] mat11 = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(a1Matrix), MatrixOperations.MatrixProduct(VectorOperations.VectorVectorTensorProduct(normalVector, dRho[0]), aMatrix));
            double[,] mat12 = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(aMatrix), MatrixOperations.MatrixProduct(VectorOperations.VectorVectorTensorProduct(dRho[0], normalVector), a1Matrix));
            double[,] mat21 = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(a1Matrix), MatrixOperations.MatrixProduct(VectorOperations.VectorVectorTensorProduct(normalVector, dRho[1]), aMatrix));
            double[,] mat22 = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(aMatrix), MatrixOperations.MatrixProduct(VectorOperations.VectorVectorTensorProduct(dRho[0], normalVector), a2Matrix));
            double[,] mat31 = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(a2Matrix), MatrixOperations.MatrixProduct(VectorOperations.VectorVectorTensorProduct(normalVector, dRho[0]), aMatrix));
            double[,] mat32 = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(aMatrix), MatrixOperations.MatrixProduct(VectorOperations.VectorVectorTensorProduct(dRho[1], normalVector), a1Matrix));
            double[,] mat41 = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(a2Matrix), MatrixOperations.MatrixProduct(VectorOperations.VectorVectorTensorProduct(normalVector, dRho[1]), aMatrix));
            double[,] mat42 = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(aMatrix), MatrixOperations.MatrixProduct(VectorOperations.VectorVectorTensorProduct(dRho[1], normalVector), a2Matrix));
            
            double[,] mat1 = MatrixOperations.MatrixAddition(mat11, mat12);
            double[,] mat2 = MatrixOperations.MatrixAddition(mat21, mat22);
            double[,] mat3= MatrixOperations.MatrixAddition(mat31, mat32);
            double[,] mat4 = MatrixOperations.MatrixAddition(mat41, mat42);

            double[,] Kr = MatrixOperations.MatrixAddition(MatrixOperations.MatrixAddition(MatrixOperations.MatrixAddition(MatrixOperations.ScalarMatrixProductNew(scalar1, mat1),
                MatrixOperations.ScalarMatrixProductNew(scalar2, mat2)),
                MatrixOperations.ScalarMatrixProductNew(scalar3, mat3)),
                MatrixOperations.ScalarMatrixProductNew(scalar4, mat4));
            return Kr;
        }
        public double[,] CreateGlobalStiffnessMatrix()
        {
            double[] ksiVector = Project(new double[2]);


            if (Math.Abs(ksiVector[0]) <= 1.05 && ksiVector[1] <= 1.05)
            {
                Tuple<double[,], double[,], double[,]> aMatrices = CalculatePositionMatrix(ksiVector[0], ksiVector[1]);
                List<double[]> dRho = SurfaceVectors(ksiVector[0], ksiVector[1]);
                double[,] m = MetricTensor(dRho);
                double[] n = NormalVector(m, dRho);
                double[] xUpdated = xUpdatedVector();
                double ksi3 = CalculatePenetration(n, aMatrices.Item1, xUpdated);

                if (ksi3 <= 0)
                {
                    double[,] Km = MainStiffnessPart(PenaltyFactor, n, aMatrices.Item1);
                    double[,] Kr = RotationalStiffnessPart(PenaltyFactor, n, aMatrices.Item1, aMatrices.Item2, aMatrices.Item3, dRho, ksi3);
                    double[,] K = MatrixOperations.MatrixAddition(Km, Kr);

                    return K;
                }
                else
                {
                    return new double[15, 15];
                }
            }
            else
            {
                return new double[15, 15];
            }
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double[] ksiVector = Project(new double[2]);
            if (Math.Abs(ksiVector[0]) <= 1.05 && ksiVector[1] <= 1.05)
            {
                Tuple<double[,], double[,], double[,]> aMatrices = CalculatePositionMatrix(ksiVector[0], ksiVector[1]);
                List<double[]> dRho = SurfaceVectors(ksiVector[0], ksiVector[1]);
                double[,] m = MetricTensor(dRho);
                double[] n = NormalVector(m, dRho);
                double[] xUpdated = xUpdatedVector();
                double ksi3 = CalculatePenetration(n, aMatrices.Item1, xUpdated);

                if (ksi3 <= 0)
                {
                    double[,] AT = MatrixOperations.Transpose(aMatrices.Item1);
                    double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
                    double[] internalGlobalForcesVector = VectorOperations.VectorScalarProductNew(AT_n, PenaltyFactor * ksi3);
                    return internalGlobalForcesVector;
                }
                else
                {
                    return new double[15];
                }
            }
            else
            {
                return new double[15];
            }
        }

        public double ClosestPointProjection()
        {
            throw new Exception("Alternative method <Project> has been used for 3D contact");
        }

        public double[,] CreateMassMatrix()
        {
            return new double[15, 15];
        }

        public double[,] CreateDampingMatrix()
        {
            return new double[15, 15];
        }
    }
}