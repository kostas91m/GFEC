using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class Triangle3 : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }

        public Triangle3(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, true, false, false, false, false };
            ElementFreedomSignature[3] = new bool[] { true, true, false, false, false, false };
            DisplacementVector = new double[6];
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
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
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
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
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
        private double[] UpdateNodalCoordinates(double[] displacementVector)
        {
            double[] updatedCoor = new double[6];
            for (int i = 1; i <= 3; i++)
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
            return finalNodes;
        }

        private Dictionary<int, double> CalculateShapeFunctions(double ksi, double ihta)
        {
            Dictionary<int, double> shapeFunctions = new Dictionary<int, double>();
            double N1 = ksi; shapeFunctions.Add(1, N1);
            double N2 = ihta; shapeFunctions.Add(2, N2);
            double N3 = 1.0 - ksi - ihta; shapeFunctions.Add(3, N3);
            return shapeFunctions;
        }

        private double[,] CalculateShapeFunctionMatrix(double ksi, double ihta)
        {
            Dictionary<int, double> shapeFunctions = CalculateShapeFunctions(ksi, ihta);
            double[,] N = new double[,]
            {
                {shapeFunctions[1], 0, shapeFunctions[2], 0, shapeFunctions[3], 0},
                {0, shapeFunctions[1], 0, shapeFunctions[2], 0, shapeFunctions[3]}
            };
            return N;
        }

        private Dictionary<string, double[]> CalculateShapeFunctionsLocalDerivatives(double[] naturalCoordinates)
        {
            double ksi = naturalCoordinates[0];
            double ihta = naturalCoordinates[1];

            double[] dN_ksi = new double[]
            {
                (1),
                (0),
                (-1),
            };

            double[] dN_ihta = new double[]
            {
                (0),
                (1),
                (-1),
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
            for (int i = 0; i < 3; i++)
            {
                jacobianMatrix[0, 0] = jacobianMatrix[0, 0] + xUpdated[k] * dN["ksi"][i];
                k = k + 2;
            }
            k = 1;
            for (int i = 0; i < 3; i++)
            {
                jacobianMatrix[0, 1] = jacobianMatrix[0, 1] + xUpdated[k] * dN["ksi"][i];
                k = k + 2;
            }

            k = 0;
            for (int i = 0; i < 3; i++)
            {
                jacobianMatrix[1, 0] = jacobianMatrix[1, 0] + xUpdated[k] * dN["ihta"][i];
                k = k + 2;
            }
            k = 1;
            for (int i = 0; i < 3; i++)
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

            for (int i = 0; i < 3; i++)
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
            double[,] Bmatrix = new double[3, 6];

            for (int i = 0; i < 3; i++)
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

        private Tuple<double[], double[]> GaussPoints(int i , int j)
        {
            double[] gaussPoints = new double[] { 0.0571041961, 0.2768430136, 0.5835904324, 0.8602401957 };
            double[] gaussPoints2 = new double[] { 0.06546699465, 0.2386486597, 0.2789904635, 0.1300560792 };
            double[] gaussWeights = new double[] { 0.1355069134, 0.2034645680, 0.1298475476, 0.0311809709 };
            double[] gaussWeights2 = new double[] { 0.1739274226, 0.3260725774, 0.3260725774, 0.1739274226 };
            double[] vectorWithPoints = new double[] { gaussPoints[i], gaussPoints2[j] };
            double[] vectorWithWeights = new double[] { gaussWeights[i], gaussWeights2[j] };
            return new Tuple<double[], double[]>(vectorWithPoints, vectorWithWeights);
        }

        public double[,] CreateGlobalStiffnessMatrix()
        {
            double[,] K = new double[6, 6];
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);

            for (int i = 0; i <= 3; i++)
            {
                for (int j = 0; j <= 3; j++)
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
            double[,] tempM = MatrixOperations.CreateDiagonalMatrix(6, 1.0);
            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            double a = Math.Pow(Math.Pow(xUpdated[0] - xUpdated[2], 2) + Math.Pow(xUpdated[1] - xUpdated[3], 2), 0.5);
            double b = Math.Pow(Math.Pow(xUpdated[0] - xUpdated[4], 2) + Math.Pow(xUpdated[1] - xUpdated[5], 2), 0.5);
            double c = Math.Pow(Math.Pow(xUpdated[2] - xUpdated[4], 2) + Math.Pow(xUpdated[3] - xUpdated[5], 2), 0.5);
            double t = (a + b + c) / 2;
            double area = Math.Pow(t * (t - a) * (t - b) * (t - c), 0.5);
            double scalar = Properties.Density * Properties.Thickness * area / 3.0;
            double[,] M = MatrixOperations.ScalarMatrixProductNew(scalar, tempM);
            double waveSpeed = Math.Sqrt(Properties.YoungMod / Properties.Density);
            return M;
        }

        public double[,] CreateDampingMatrix()
        {
            return new double[6, 6];
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double[] F = new double[6];
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio); 

            for (int i = 0; i <= 3; i++)
            {
                for (int j = 0; j <= 3; j++)
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

