using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class ContactNtS2DTh : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        private double PenaltyFactor { get; set; }
        private double ContactArea { get; set; }
        private double ContactPressure { get; set; }
        private double ContactThermalConductivity { get; set; }
        private double SurfaceRoughness { get; set; }
        private double YieldStrength { get; set; }
        private double Dx1 { get; set; }
        private double Dx2 { get; set; }
        private double Dx { get; set; }

        public ContactNtS2DTh(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, false, false, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, false, false, false, false, false };
            ElementFreedomSignature[3] = new bool[] { true, false, false, false, false, false };
            DisplacementVector = new double[3];
            ContactArea = properties.SectionArea;
            ContactPressure = properties.ContactForceValue / properties.SectionArea;
            SurfaceRoughness = properties.SurfaceRoughness;
            ContactThermalConductivity = properties.ContactThermalConductivity;
            YieldStrength = properties.YieldStrength;
            Dx = properties.Dx;
            Dx1 = properties.Dx1;
            Dx2 = properties.Dx2;
        }

        public double ClosestPointProjection()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }

        public Dictionary<int, INode> NodesAtFinalState()
        {
            throw new Exception("Method not implemenented");
        }

        public double CalculateConductivity()
        {
            //double k = 19.2;
            double m = 0.1259 * Math.Pow(SurfaceRoughness * Math.Pow(10, 6), 0.402);
            //double c1 = 6271.0 * Math.Pow(10, 6);//SS
            //double c2 = -0.229;//SS
            double c1 = 1024691.406 * Math.Pow(10, 3);//154237.1846;//180460.6748 * Math.Pow(10, 3);//Yovanovitch5//200577.997 * Math.Pow(10, 3);//Yovanovitch4//236818.1659 * Math.Pow(10, 3);//Yovanovitch3//267224.3367 * Math.Pow(10, 3);//Yovanovitch2//453323.5643 * Math.Pow(10, 3);//Yovanovitch1
            double c2 = -0.745;//-0.118//Yovanovitch6// - 0.17;//Yovanovitch5//-0.205;//Yovanovitch4// -0.26;//Yovanovitch3// -0.3;//Yovanovitch2//-0.475;//Yovanovitch1
            //double sigma = 0.478 * Math.Pow(10, -6);
            double cc = (1.25 * ContactThermalConductivity * m / SurfaceRoughness) * Math.Pow((ContactPressure / c1) * Math.Pow(1.6177 * 1000000 * SurfaceRoughness / m, -c2), 0.95 / (1 + 0.0711 * c2));
            //double cc = 1.25 * ContactThermalConductivity * (m / (SurfaceRoughness)) * Math.Pow(ContactPressure / (3.0 * YieldStrength), 0.95);// * 2.582260191 * Math.Pow(10, -4);
            double cH = cc * ContactArea;
            return cH;
        }

        public double CalculateCoefficient1()
        {
            double W1;
            W1 = 1.0 - Dx1;  //(Dx - Dx1) / Dx;
            return W1;
        }
        public double CalculateCoefficient2()
        {
            double W2;
            double W1 = CalculateCoefficient1();
            W2 = 1.0 - W1;
            //W2 = 1.0 - Dx2 / Dx;
            return W2;
        }
        private double CalculateTemperatureJump()
        {
            double W1 = CalculateCoefficient1();
            double W2 = CalculateCoefficient2();
            double theta1 = 0.0;
            double theta2 = 0.0;
            double theta3 = 0.0;
            double gH = (theta3 + DisplacementVector[2]) - (W1 * (theta1 + DisplacementVector[0]) + W2 * (theta2 + DisplacementVector[1]));
            return gH;
        }
        private double CalculateTemperatureJump1()
        {
            //double W1 = CalculateCoefficient1();
            double theta1 = 0.0;
            //double theta2 = 0.0;
            double theta3 = 0.0;
            double gH1 = (theta3 + DisplacementVector[2]) - (theta1 + DisplacementVector[0]);
            return gH1;
        }
        private double CalculateTemperatureJump2()
        {
            //double W2 = CalculateCoefficient2();
            //double theta1 = 0.0;
            double theta2 = 0.0;
            double theta3 = 0.0;
            double gH2 = (theta3 + DisplacementVector[2]) - (theta2 + DisplacementVector[1]);
            return gH2;
        }
        public double[,] CreateGlobalStiffnessMatrix()
        {
            double cH = CalculateConductivity();
            double W1 = CalculateCoefficient1();
            double W2 = CalculateCoefficient2();
            double[,] stiffMatrix = new double[,]
            {
                { W1 * cH, 0.0, -W1 * cH},
                { 0.0, W2 * cH, -W2 * cH},
                { -W1 * cH, -W2 * cH, 1.0 * cH},
            };
            return stiffMatrix;
            //{
            //    {0.50 * W1 * cH, 0.50 * W2 * cH, -0.50 * cH},
            //    {0.50 * W1 * cH, 0.50 * W2 * cH, -0.50 * cH},
            //    {-W1 * cH, -W2 * cH, 1.0 * cH},
            //};
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double cH = CalculateConductivity();
            double gH = CalculateTemperatureJump();
            double gH1 = CalculateTemperatureJump1();
            double gH2 = CalculateTemperatureJump2();
            double W1 = CalculateCoefficient1();
            double W2 = CalculateCoefficient2();
            double[] heatFlux = new double[]
            {
                -W1 * cH * gH1 ,
                -W2 * cH * gH2 ,
                1.0 * cH * gH
            };
            return heatFlux;
        }

        public double[,] CreateMassMatrix()
        {
            return new double[6, 6];
        }

        public double[,] CreateDampingMatrix()
        {
            return new double[6, 6];
        }
    }
}

