using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class SecondOrderContactNtS2D : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        private double PenaltyFactor { get; set; }

        public SecondOrderContactNtS2D(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[3] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[4] = new bool[] { true, true, false, false, false, false };
            DisplacementVector = new double[8];
            PenaltyFactor = properties.YoungMod * 10.0;// εN / Ε is to be added to element properties
        }
        public Dictionary<int, INode> NodesAtFinalState()
        {
            Dictionary<int, INode> finalNodes = new Dictionary<int, INode>();
            finalNodes[1] = new Node(Nodes[1].XCoordinate + DisplacementVector[0], Nodes[1].YCoordinate + DisplacementVector[1]);
            finalNodes[2] = new Node(Nodes[2].XCoordinate + DisplacementVector[2], Nodes[2].YCoordinate + DisplacementVector[3]);
            finalNodes[3] = new Node(Nodes[3].XCoordinate + DisplacementVector[4], Nodes[3].YCoordinate + DisplacementVector[5]);
            finalNodes[4] = new Node(Nodes[4].XCoordinate + DisplacementVector[6], Nodes[4].YCoordinate + DisplacementVector[7]);
            return finalNodes;
        }
        private double[] NodalXUpdated()
        {
            double[] x = new double[] { Nodes[1].XCoordinate + DisplacementVector[0],
            Nodes[1].YCoordinate + DisplacementVector[1],
            Nodes[2].XCoordinate + DisplacementVector[2],
            Nodes[2].YCoordinate + DisplacementVector[3],
            Nodes[3].XCoordinate + DisplacementVector[4],
            Nodes[3].YCoordinate + DisplacementVector[5],
            Nodes[4].XCoordinate + DisplacementVector[6],
            Nodes[4].YCoordinate + DisplacementVector[7],
            };
            return x;
        }
        public List<double[]> GetStressVector()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0, 0.0 });
            return l;
        }
        public List<double[]> GetStrainVector()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0, 0.0 });
            return l;
        }

        public List<double[]> GetGaussPointsInPhysicalSpace()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0 });
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
        public double ClosestPointProjection()
        {
              throw new Exception("Alternative method <Project> has been used for higher order elements");
        }

        private Tuple<double[,], double[,], double[,]> CalculatePositionMatrix(double ksi1)
        {
            double N1 = 1.0 / 2.0 * (Math.Pow(ksi1, 2.0) - ksi1);
            double N2 = 1.0 - Math.Pow(ksi1, 2.0);
            double N3 = 1.0 / 2.0 * (Math.Pow(ksi1, 2.0) + ksi1);
            double dN1 = ksi1 - 0.5;
            double dN2 = -2.0 * ksi1;
            double dN3 = ksi1 + 0.5;
            double dN21 = 1.0;
            double dN22 = -2.0;
            double dN23 = 1.0;
            double[,] aMatrix = new double[,]
                {
                    { -N1 ,0.0 ,-N2 ,0.0 ,-N3 , 0.0 , 1.0 ,0.0 },
                    {0.0, -N1 , 0.0 ,-N2, 0.0, -N3, 0.0, 1.0 }
                };

            double[,] daMatrix = new double[,]
                {
                    { -dN1 ,0.0 ,-dN2 ,0.0 ,-dN3, 0.0, 0.0, 0.0 },
                    {0.0, -dN1 , 0.0 ,-dN2, 0.0, -dN3, 0.0, 0.0 }
                };
            double[,] da2Matrix = new double[,]
    {
                    { -dN21 ,0.0 ,-dN22 ,0.0 ,-dN23, 0.0, 0.0, 0.0 },
                    {0.0, -dN21 , 0.0 ,-dN22, 0.0, -dN23, 0.0, 0.0 }
    };
            return new Tuple<double[,], double[,], double[,]>(aMatrix, daMatrix, da2Matrix);
        }

        private Tuple<double[], double, double[], double[], double> MasterSegmentGeometry(double[,] daMatrix, double[,] da2Matrix)
        {
            double Xm1 = Nodes[1].XCoordinate + DisplacementVector[0];
            double Ym1 = Nodes[1].YCoordinate + DisplacementVector[1];
            double Xm2 = Nodes[2].XCoordinate + DisplacementVector[2];
            double Ym2 = Nodes[2].YCoordinate + DisplacementVector[3];
            double Xm3 = Nodes[3].XCoordinate + DisplacementVector[4];
            double Ym3 = Nodes[3].YCoordinate + DisplacementVector[5];
            double Xs = Nodes[4].XCoordinate + DisplacementVector[6];
            double Ys = Nodes[4].YCoordinate + DisplacementVector[7];

            double[] xupd = new double[] { -Xm1, -Ym1, -Xm2, -Ym2, -Xm3, -Ym3, - Xs, -Ys };
            double[] surfaceVector = VectorOperations.MatrixVectorProduct(daMatrix, xupd);
            double[] surfaceVectorDerivative = VectorOperations.MatrixVectorProduct(da2Matrix, xupd);

            double detm = VectorOperations.VectorDotProduct(surfaceVector, surfaceVector);
            double m11 = 1.0 / detm;
            double[] vector = new double[] { -surfaceVector[1], surfaceVector[0]};
            double scalarCoef = 1.0 / (Math.Sqrt(detm));
            double[] normalUnitVec = VectorOperations.VectorScalarProductNew(vector, scalarCoef);
            double[] tangentVector = VectorOperations.VectorScalarProductNew(surfaceVector, scalarCoef);
            double scalarCoef2 = Math.Pow(m11, 2.0);
            double curvatureTensor = scalarCoef2 * VectorOperations.VectorDotProduct(surfaceVectorDerivative, normalUnitVec);

            return new Tuple<double[], double, double[], double[], double>(surfaceVector, m11, normalUnitVec, tangentVector, curvatureTensor);
        }

        private double CalculatePenetration(double[,] aMatrix, double[] n)
        {
            double[,] AT = MatrixOperations.Transpose(aMatrix);
            double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
            double[] xupd = new double[] {
                Nodes[1].XCoordinate + DisplacementVector[0],
                Nodes[1].YCoordinate + DisplacementVector[1],
                Nodes[2].XCoordinate + DisplacementVector[2],
                Nodes[2].YCoordinate + DisplacementVector[3],
                Nodes[3].XCoordinate + DisplacementVector[4],
                Nodes[3].YCoordinate + DisplacementVector[5],
                Nodes[4].XCoordinate + DisplacementVector[6],
                Nodes[4].YCoordinate + DisplacementVector[7]
            };
            double normalGap = VectorOperations.VectorDotProduct(xupd, AT_n);
            return normalGap;
        }
        private double CalculateDeltaKsi(double[] masterSlaveRelativeVector, double[] surfaceVector, double[] surfaceVectorDerivative)
        {
            double scalar1 = VectorOperations.VectorDotProduct(surfaceVector, masterSlaveRelativeVector);
            double scalar2 = VectorOperations.VectorDotProduct(surfaceVectorDerivative, masterSlaveRelativeVector) -
                             VectorOperations.VectorDotProduct(surfaceVector, surfaceVector);

            double deltaKsi = -scalar1 / scalar2;
            return deltaKsi;
        }
        private double Project(double ksi1Initial)
        {
            int maxIterations = 1000;
            double tol = Math.Pow(10.0, -6.0);
            double deltaKsi = 0.0;
            double ksi = ksi1Initial;
            double[] xUpdated = NodalXUpdated();
            for (int i = 1; i <= maxIterations; i++)
            {
                Tuple<double[,], double[,], double[,]> aMatrices = CalculatePositionMatrix(ksi);
                //double[] slavePositionVector = new double[] { xUpdated[6], xUpdated[7] };
                double[] masterSlaveRelativeVector = VectorOperations.MatrixVectorProduct(aMatrices.Item1, xUpdated);
                double[] surfaceVector = VectorOperations.VectorScalarProduct(VectorOperations.MatrixVectorProduct(aMatrices.Item2, xUpdated), -1);
                double[] surfaceVectorDerivative = VectorOperations.VectorScalarProduct(VectorOperations.MatrixVectorProduct(aMatrices.Item3, xUpdated), -1);

                deltaKsi = CalculateDeltaKsi(masterSlaveRelativeVector, surfaceVector, surfaceVectorDerivative);
                ksi += deltaKsi;
                if (Math.Abs(deltaKsi) <= tol)
                {
                    break;
                }
            }
            if (Math.Abs(deltaKsi) > tol)
            {
                throw new Exception("CPP not found in current iterations");
            }
            else
            {
                return ksi;

            }
        }

        private double[,] CalculateMainStiffnessPart(double ksi1, double[] n)
        {
            double[,] mainStiffnessMatrix = new double[8, 8];
            Tuple<double[,], double[,], double[,]> positionMatrices = CalculatePositionMatrix(ksi1);
            double[,] A = positionMatrices.Item1;
            double[,] nxn = VectorOperations.VectorVectorTensorProduct(n, n);
            double[,] nxn_A = MatrixOperations.MatrixProduct(nxn, A);
            double[,] AT_nxn_A = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(A), nxn_A);
            mainStiffnessMatrix = MatrixOperations.ScalarMatrixProductNew(PenaltyFactor, AT_nxn_A);
            return mainStiffnessMatrix;
        }

        private double[,] CalculateRotationalStiffnessPart(double[,] A, double[,] dA, double[] n, double ksi3, double m11, double[] dRho)
        {
            double coef = PenaltyFactor * ksi3 * m11;
            double[,] rotationalPart;
            double[,] n_x_dRho = VectorOperations.VectorVectorTensorProduct(n, dRho);
            double[,] dRho_x_n = VectorOperations.VectorVectorTensorProduct(dRho, n);
            double[,] firstTerm = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(dA),
                                                                    MatrixOperations.MatrixProduct(n_x_dRho, A)
                                                                    );
            double[,] secondTerm = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(A),
                                                                    MatrixOperations.MatrixProduct(dRho_x_n, dA)
                                                                    );
            rotationalPart = MatrixOperations.ScalarMatrixProductNew(
                                                                        coef,
                                                                        MatrixOperations.MatrixAddition(firstTerm, secondTerm)
                                                                        );
            return rotationalPart;
        }
        private double[,] CalculateCurvatureStiffnessPart(double[,] A, double ksi3, double m11, double[] dRho, double h11)
        {
            double coef = PenaltyFactor * ksi3 * m11 * h11;
            double[,] curvaturePart;
            double[,] dRho_x_dRho = VectorOperations.VectorVectorTensorProduct(dRho, dRho);
            double[,] Matrix = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(A),
                                                                    MatrixOperations.MatrixProduct(dRho_x_dRho, A)
                                                                    );
            curvaturePart = MatrixOperations.ScalarMatrixProductNew(
                                                                        coef,
                                                                        Matrix
                                                                        );
            return curvaturePart;
        }
        public double[,] CreateGlobalStiffnessMatrix()
        {
            double ksi1 = Project(0.0);
            if (Math.Abs(ksi1) <= 1.05)
            {
                Tuple<double[,], double[,], double[,]> positionMatrices = CalculatePositionMatrix(ksi1);
                double[,] aMatrix = positionMatrices.Item1;
                double[,] daMatrix = positionMatrices.Item2;
                double[,] da2Matrix = positionMatrices.Item3;

                Tuple<double[], double, double[], double[], double> surfaceCharacteristics = MasterSegmentGeometry(daMatrix, da2Matrix);
                double m11 = surfaceCharacteristics.Item2;
                double[] dRho = surfaceCharacteristics.Item1;
                double[] n = surfaceCharacteristics.Item3;
                double h11 = surfaceCharacteristics.Item5;
                double ksi3 = CalculatePenetration(aMatrix, n);
                if (ksi3 <= 0)
                {
                    double[,] mainPart = CalculateMainStiffnessPart(ksi1, n);
                    double[,] rotationalPart = CalculateRotationalStiffnessPart(aMatrix, daMatrix, n, ksi3, m11, dRho);
                    double[,] curvaturePart = CalculateCurvatureStiffnessPart(aMatrix, ksi3, m11, dRho, h11);
                    double[,] globalStiffnessMatrix = MatrixOperations.MatrixAddition(MatrixOperations.MatrixAddition(mainPart, rotationalPart), curvaturePart);
                    return globalStiffnessMatrix;
                }
                else
                {
                    double[,] globalStifnessMatrix = new double[8, 8];
                    return globalStifnessMatrix;
                }
            }
            else
            {
                double[,] globalStifnessMatrix = new double[8, 8];
                return globalStifnessMatrix;
            }
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double ksi1 = Project(0.0);
            if (Math.Abs(ksi1) <= 1.05)
            {
                Tuple<double[,], double[,], double[,]> positionMatrices = CalculatePositionMatrix(ksi1);
                double[,] aMatrix = positionMatrices.Item1;
                double[,] daMatrix = positionMatrices.Item2;
                double[,] da2Matrix = positionMatrices.Item3;

                Tuple<double[], double, double[], double[], double> surfaceCharacteristics = MasterSegmentGeometry(daMatrix, da2Matrix);
                double[] n = surfaceCharacteristics.Item3;
                double ksi3 = CalculatePenetration(aMatrix, n);
                if (ksi3 <= 0)
                {
                    double[,] AT = MatrixOperations.Transpose(aMatrix);
                    double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
                    double[] internalGlobalForcesVector = VectorOperations.VectorScalarProductNew(AT_n, PenaltyFactor * ksi3);
                    return internalGlobalForcesVector;
                }
                else
                {
                    double[] internalGlobalForcesVector = new double[8];
                    return internalGlobalForcesVector;
                }
            }
            else
            {
                double[] internalGlobalForcesVector = new double[8];
                return internalGlobalForcesVector;
            }
        }

        public double[,] CreateMassMatrix()
        {
            return new double[8, 8];
        }

        public double[,] CreateDampingMatrix()
        {
            return new double[8, 8];
        }
    }
}

