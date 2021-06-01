using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class Quad8 : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        public double poisson { get; set; }
        //private double thickness = 1.0; //To be included in Element Properties
        //private double density = 1.0; //To be included in Element Properties

        public Quad8(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[3] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[4] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[5] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[6] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[7] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[8] = new bool[] { true, true, false, false, false, false };
            DisplacementVector = new double[16];
        }

        public double ClosestPointProjection()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public List<double[]> GetStressVector()
        {
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            List<double[]> stressVectors = new List<double[]>();
            List<double[]> strainVectors = GetStrainVector();
            foreach (var v in strainVectors)
            {
                double[] stressV = CalculateStressVector(E, v);
                stressVectors.Add(stressV);
            }
            return stressVectors;
        }
        public List<double[]> GetStrainVector()
        {

            List<double[]> strainVectors = new List<double[]>();
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    double[] gP = GaussPoints(i, j).Item1;
                    Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                    double[,] J = CalculateJacobian(localdN);
                    double[,] invJ = CalculateInverseJacobian(J).Item1;
                    double detJ = CalculateInverseJacobian(J).Item2;
                    Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                    double[,] B = CalculateBMatrix(globaldN);
                    double[] strainV = CalculateStrainsVector(B);
                    strainVectors.Add(strainV);
                }
            }
            return strainVectors;
        }
        public List<double[]> GetGaussPointsInPhysicalSpace()
        {
            List<double[]> gaussPoints = new List<double[]>();
            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    double[] gP = GaussPoints(i, j).Item1;
                    double[] gaussPoint = VectorOperations.MatrixVectorProduct(CalculateShapeFunctionMatrix(gP[0], gP[1]), xUpdated);
                    gaussPoints.Add(gaussPoint);
                }
            }
            return gaussPoints;
        }
        public List<double[]> GetStressFromElementsNodes()
        {
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            List<double[]> stressVectors = new List<double[]>();
            List<double[]> strainVectors = GetStrainFromElementsNodes();
            foreach (var v in strainVectors)
            {
                double[] stressV = CalculateStressVector(E, v);
                stressVectors.Add(stressV);
            }
            return stressVectors;
        }
        public List<double[]> GetStrainFromElementsNodes()
        {
            List<double[]> strainVectors = new List<double[]>();
            double[] ksi = new double[] { -1.0, 0.0, 1.0, 1.0, 1.0, 0.0, -1.0, - 1.0 };
            double[] ihta = new double[] { -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 1.0, 0.0 };
            for (int i = 0; i < 8; i ++)
            {
                double[] node = new double[] { ksi[i], ihta[i] };
                Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(node);
                double[,] J = CalculateJacobian(localdN);
                double[,] invJ = CalculateInverseJacobian(J).Item1;
                double detJ = CalculateInverseJacobian(J).Item2;
                Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                double[,] B = CalculateBMatrix(globaldN);
                double[] strainV = CalculateStrainsVector(B);
                strainVectors.Add(strainV);
            }
            return strainVectors;
        }
        private double[] UpdateNodalCoordinates(double[] displacementVector)
        {
            double[] updatedCoor = new double[16];
            for (int i = 1; i <= 8; i++)
            {
                updatedCoor[2 * i - 2] = Nodes[i].XCoordinate + displacementVector[2 * i - 2];
                updatedCoor[2 * i - 1] = Nodes[i].YCoordinate + displacementVector[2 * i - 1];
            }
            return updatedCoor;
        }

        public Dictionary<int, INode> NodesAtFinalState()
        {
            Dictionary<int, INode> finalNodes = new Dictionary<int, INode>();
            finalNodes[1] = new Node(Nodes[1].XCoordinate + DisplacementVector[0], Nodes[1].YCoordinate + DisplacementVector[1]);
            finalNodes[2] = new Node(Nodes[2].XCoordinate + DisplacementVector[2], Nodes[2].YCoordinate + DisplacementVector[3]);
            finalNodes[3] = new Node(Nodes[3].XCoordinate + DisplacementVector[4], Nodes[3].YCoordinate + DisplacementVector[5]);
            finalNodes[4] = new Node(Nodes[4].XCoordinate + DisplacementVector[6], Nodes[4].YCoordinate + DisplacementVector[7]);
            finalNodes[5] = new Node(Nodes[5].XCoordinate + DisplacementVector[8], Nodes[5].YCoordinate + DisplacementVector[9]);
            finalNodes[6] = new Node(Nodes[6].XCoordinate + DisplacementVector[10], Nodes[6].YCoordinate + DisplacementVector[11]);
            finalNodes[7] = new Node(Nodes[7].XCoordinate + DisplacementVector[12], Nodes[7].YCoordinate + DisplacementVector[13]);
            finalNodes[8] = new Node(Nodes[8].XCoordinate + DisplacementVector[14], Nodes[8].YCoordinate + DisplacementVector[15]);
            return finalNodes;
        }

        private Dictionary<int, double> CalculateShapeFunctions(double ksi, double ihta)
        {
            Dictionary<int, double> shapeFunctions = new Dictionary<int, double>();
            double N1 = -1.0 / 4.0 * (1 - ksi) * (1 - ihta) * (1 + ksi + ihta); shapeFunctions.Add(1, N1);
            double N2 = 1.0 / 2.0 * (1 + ksi) * (1 - ksi) * (1 - ihta); shapeFunctions.Add(2,N2);
            double N3 = -1.0 / 4.0 * (1 + ksi) * (1 - ihta) * (1 - ksi + ihta); shapeFunctions.Add(3, N3);
            double N4 = 1.0 / 2.0 * (1 + ksi) * (1 + ihta) * (1 - ihta); shapeFunctions.Add(4, N4);
            double N5 = -1.0 / 4.0 * (1 + ksi) * (1 + ihta) * (1 - ksi - ihta); shapeFunctions.Add(5, N5);
            double N6 = 1.0 / 2.0 * (1 - ksi) * (1 + ksi) * (1 + ihta); shapeFunctions.Add(6, N6);
            double N7 = -1.0 / 4.0 * (1 - ksi) * (1 + ihta) * (1 + ksi - ihta); shapeFunctions.Add(7, N7);
            double N8 = 1.0 / 2.0 * (1 - ksi) * (1 + ihta) * (1 - ihta); shapeFunctions.Add(8, N8);
            return shapeFunctions;
        }

        private double[,] CalculateShapeFunctionMatrix(double ksi, double ihta)
        {
            Dictionary<int, double> shapeFunctions = CalculateShapeFunctions(ksi, ihta);
            double[,] N = new double[,]
            {
                {shapeFunctions[1], 0, shapeFunctions[2], 0, shapeFunctions[3], 0, shapeFunctions[4], 0, shapeFunctions[5], 0, shapeFunctions[6], 0, shapeFunctions[7], 0, shapeFunctions[8], 0  },
                {0, shapeFunctions[1], 0, shapeFunctions[2], 0, shapeFunctions[3], 0, shapeFunctions[4], 0, shapeFunctions[5], 0, shapeFunctions[6], 0, shapeFunctions[7], 0, shapeFunctions[8] }
            };
            return N;
        }

        private Dictionary<string, double[]> CalculateShapeFunctionsLocalDerivatives(double[] naturalCoordinates)
        {
            double ksi = naturalCoordinates[0];
            double ihta = naturalCoordinates[1];

            double[] dN_ksi = new double[]
            {
                (1.0/4.0*(1-ihta)*(2*ksi + ihta)),
                ((ihta-1)*ksi),
                (1.0/4.0*(ihta - 1)*(ihta - 2*ksi )),
                (-1.0/2.0*(ihta - 1) * (1 + ihta)),
                (1.0/4.0*(1 + ihta)*(2*ksi + ihta)),
                (-(ihta+1)*ksi),
                (-1.0/4.0*(1+ihta)*(-2*ksi + ihta)),
                (1.0/2.0*(1 + ihta) * (ihta - 1)),
            };
            double[] dN_ihta = new double[]
            {
                (1.0/4.0*(1-ksi)*(ksi + 2*ihta)),
                (1.0/2.0*(Math.Pow(ksi,2) - 1)),
                (-1.0/4.0*(1+ksi)*(ksi - 2*ihta)),
                (-(ksi+1)*ihta),
                (1.0/4.0*(1+ksi)*(ksi + 2*ihta)),
                (-1.0/2.0*(1 + ksi) * (ksi - 1)),
                (1.0/4.0*(1-ksi)*(-ksi + 2*ihta)),
                ((-1 + ksi)*ihta),
            };
            Dictionary<string, double[]> dN = new Dictionary<string, double[]>();
            dN.Add("ksi", dN_ksi);
            dN.Add("ihta", dN_ihta);
            return dN;
        }

        private double[,] CalculateJacobian(Dictionary<string, double[]> dN)
        {
            double[,] jacobianMatrix = new double[2, 2];

            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);

            int k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[0, 0] = jacobianMatrix[0, 0] + xUpdated[k] * dN["ksi"][i];
                k = k + 2;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[0, 1] = jacobianMatrix[0, 1] + xUpdated[k] * dN["ksi"][i];
                k = k + 2;
            }

            k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[1, 0] = jacobianMatrix[1, 0] + xUpdated[k] * dN["ihta"][i];
                k = k + 2;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[1, 1] = jacobianMatrix[1, 1] + xUpdated[k] * dN["ihta"][i];
                k = k + 2;
            }

            return jacobianMatrix;
        }

        private Tuple<double[,], double> CalculateInverseJacobian(double[,] jacobianMatrix)
        {
            double[,] jacobianInverseMatrix = new double[2, 2];

            double detj = jacobianMatrix[0, 0] * jacobianMatrix[1, 1] - jacobianMatrix[0, 1] * jacobianMatrix[1, 0];

            jacobianInverseMatrix[0, 0] = jacobianMatrix[1, 1] / detj;
            jacobianInverseMatrix[0, 1] = -jacobianMatrix[0, 1] / detj;
            jacobianInverseMatrix[1, 0] = -jacobianMatrix[1, 0] / detj;
            jacobianInverseMatrix[1, 1] = jacobianMatrix[0, 0] / detj;

            return new Tuple<double[,], double>(jacobianInverseMatrix, detj);
        }

        private Dictionary<int, double[]> CalculateShapeFunctionsGlobalDerivatives(Dictionary<string, double[]> dN, double[,] Jinv)
        {
            Dictionary<int, double[]> dNg = new Dictionary<int, double[]>();

            for (int i = 0; i < 8; i++)
            {
                double[] dNlocal = new double[] { dN["ksi"][i], dN["ihta"][i] };
                double[] dNglobal = VectorOperations.MatrixVectorProduct(Jinv, dNlocal);
                dNg.Add(i, dNglobal);
            }
            return dNg;
        }

        private double[] CalculateStrainsVector(double[,] Bmatrix)
        {
            double[] strains = VectorOperations.MatrixVectorProduct(Bmatrix, DisplacementVector);
            return strains;
        }

        private double[,] CalculateBMatrix(Dictionary<int, double[]> dNglobal)
        {
            double[,] Bmatrix = new double[3, 16];

            for (int i = 0; i < 8; i++)
            {
                Bmatrix[0, i * 2] = dNglobal[i][0];
                Bmatrix[1, i * 2 + 1] = dNglobal[i][1];
                Bmatrix[2, i * 2] = dNglobal[i][1];
                Bmatrix[2, i * 2 + 1] = dNglobal[i][0];
            }
            return Bmatrix;
        }

        private double[,] CalculateStressStrainMatrix(double E, double v)
        {
            double[,] Ematrix = new double[3, 3];
            //v = 0.30;
            //v = Properties.PoissonRatio;
            double Ehat = E / ((1.0 - Math.Pow(v, 2)));

            Ematrix[0, 0] = Ehat;
            Ematrix[0, 1] = Ehat * v;
            Ematrix[1, 0] = Ehat * v;
            Ematrix[1, 1] = Ehat;
            Ematrix[2, 2] = Ehat * (1.0 / 2.0) * (1.0 - v);

            return Ematrix;
        }

        private double[] CalculateStressVector(double[,] E, double[] strain)
        {
            double[] stressVector = VectorOperations.MatrixVectorProduct(E, strain);
            return stressVector;
        }

        private Tuple<double[], double[]> GaussPoints(int i, int j)
        {
            //double[] gaussPoints = new double[] { -1.0 / Math.Sqrt(3), 1.0 / Math.Sqrt(3) };
            //double[] gaussWeights = new double[] { 1.0, 1.0 };
            double[] gaussPoints = new double[] { -0.77459, 0.0, 0.77459 };
            double[] gaussWeights = new double[] { 0.55555, 0.88888, 0.55555 };
            double[] vectorWithPoints = new double[] { gaussPoints[i], gaussPoints[j] };
            double[] vectorWithWeights = new double[] { gaussWeights[i], gaussWeights[j] };
            return new Tuple<double[], double[]>(vectorWithPoints, vectorWithWeights);
        }

        public double[,] CreateGlobalStiffnessMatrix()
        {
            double[,] K = new double[16, 16];
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    double[] gP = GaussPoints(i, j).Item1;
                    double[] gW = GaussPoints(i, j).Item2;
                    Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                    double[,] J = CalculateJacobian(localdN);
                    double[,] invJ = CalculateInverseJacobian(J).Item1;
                    double detJ = CalculateInverseJacobian(J).Item2;
                    Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                    double[,] B = CalculateBMatrix(globaldN);
                    K = MatrixOperations.MatrixAddition(K, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * Properties.Thickness,
                        MatrixOperations.MatrixProduct(MatrixOperations.Transpose(B), MatrixOperations.MatrixProduct(E, B))));
                }
            }
            return K;
        }

        public double[,] CreateMassMatrix()
        {
            double[,] M = new double[16, 16];
            double Mtot = new double();
            //double scalar = 28.24;
            //double[,] M = MatrixOperations.CreateDiagonalMatrix(8, scalar);

            double[,] consinstentMass = new double[16, 16];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    double[] gP = GaussPoints(i, j).Item1;
                    double[] gW = GaussPoints(i, j).Item2;
                    Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                    double[,] J = CalculateJacobian(localdN);
                    double detJ = CalculateInverseJacobian(J).Item2;
                    double[,] Nmatrix = CalculateShapeFunctionMatrix(gP[0], gP[1]);
                    consinstentMass = MatrixOperations.MatrixAddition(consinstentMass, MatrixOperations.ScalarMatrixProductNew(Properties.Density * Properties.Thickness * detJ * gW[0] * gW[1],
                        MatrixOperations.MatrixProduct(MatrixOperations.Transpose(Nmatrix), Nmatrix)));
                    Mtot += Properties.Density * Properties.Thickness * detJ * gW[0] * gW[1];
                }
            }
            double c = Mtot / MatrixOperations.Trace(consinstentMass);
            //M = consinstentMass;
            for (int i = 0; i <= 15; i++)
            {
                M[i, i] = c * consinstentMass[i, i];
            }
            //for (int i = 0; i <= 15; i++)
            //{
            //    for (int j = 0; j <= 15; j++)
            //    {
            //        M[i, i] += Math.Abs(consinstentMass[i, j]);
            //    }
            //}
            //-------------------------------------------------------------------
            //double[,] tempM = MatrixOperations.CreateDiagonalMatrix(8, 1.0);
            //double length = 0.3;
            //double scalar = Properties.Density * Properties.Thickness * length * (length / 3.0) / 4.0;
            //double[,] M = MatrixOperations.ScalarMatrixProductNew(scalar, tempM);

            //double waveSpeed = Math.Sqrt(Properties.YoungMod / Properties.Density);
            //double deltatCritical = length * Math.Sqrt(1.0 - 0.33) / waveSpeed;


            //--------------------------------------------------------------
            //for (int i = 0; i < 2; i++)
            //{
            //    for (int j = 0; j < 2; j++)
            //    {
            //        double[] gP = GaussPoints(i, j).Item1;
            //        double[] gW = GaussPoints(i, j).Item2;
            //        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
            //        double[,] J = CalculateJacobian(localdN);
            //        double[,] invJ = CalculateInverseJacobian(J).Item1;
            //        double detJ = CalculateInverseJacobian(J).Item2;
            //        double[,] Nmatrix = CalculateShapeFunctionMatrix(gP[i], gP[j]);
            //        M = MatrixOperations.MatrixAddition(M, MatrixOperations.ScalarMatrixProductNew(Properties.Density * Properties.Thickness * detJ * gW[i] * gW[j],
            //            MatrixOperations.MatrixProduct(MatrixOperations.Transpose(Nmatrix), Nmatrix)));
            //    }
            //}

            //--------------------------------------------------------

            //for (int i = 0; i < 8; i++)
            //{
            //    M[i, i] = 4.0;
            //}

            //for (int i = 0; i < 6; i++)
            //{
            //    M[i, i + 2] = 2.0;
            //    M[i + 2, i] = 2.0;
            //}

            //for (int i = 0; i < 4; i++)
            //{
            //    M[i, i + 4] = 1.0;
            //    M[i + 4, i] = 1.0;
            //}

            //for (int i = 0; i < 2; i++)
            //{
            //    M[i, i + 6] = 2.0;
            //    M[i + 6, i] = 2.0;
            //}

            //M = MatrixOperations.ScalarMatrixProductNew(0.67 * 0.8 * Properties.Density * Properties.Thickness / 32, M);
            //MatrixOperations.PrintMatrix(M);

            return M;
        }

        public double[,] CreateDampingMatrix()
        {
            return new double[16, 16];
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double[] F = new double[16];
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio); //needs fixing in poisson v

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    double[] gP = GaussPoints(i, j).Item1;
                    double[] gW = GaussPoints(i, j).Item2;
                    Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                    double[,] J = CalculateJacobian(localdN);
                    double[,] invJ = CalculateInverseJacobian(J).Item1;
                    double detJ = CalculateInverseJacobian(J).Item2;
                    Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                    double[,] B = CalculateBMatrix(globaldN);
                    double[] strainVector = CalculateStrainsVector(B);
                    double[] stressVector = CalculateStressVector(E, strainVector);
                    F = VectorOperations.VectorVectorAddition(F, VectorOperations.VectorScalarProductNew(
                        VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(B), stressVector), detJ * gW[0] * gW[1] * Properties.Thickness));
                }
            }
            return F;
        }
    }
}

