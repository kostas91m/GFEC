using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    class ContactNtN2DTh : IElement
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

        public ContactNtN2DTh(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, false, false, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, false, false, false, false, false };
            DisplacementVector = new double[2];
            ContactArea = properties.SectionArea;
            ContactPressure = properties.ContactForceValue / properties.SectionArea;
        }

        public Dictionary<int, INode> NodesAtFinalState()
        {
            throw new Exception("Method not implemenented");
        }

        private double CalculateConductivity()
        {
            double cc = 1.25 * Math.Pow(ContactPressure / (3.0 * 250.0 * Math.Pow(10, 6)), 0.95);//19.2;
            double cH = cc * ContactArea;
            return cH;
        }

        private double CalculateTemperatureJump()
        {
            double theta1 = 0.0;
            double theta2 = 0.0;
            double gH = (theta2 + DisplacementVector[1]) - (theta1 + DisplacementVector[0]);
            return gH;
        }

        public double[,] CreateGlobalStiffnessMatrix()
        {
            double cH = CalculateConductivity();
            double[,] stiffMatrix = new double[,]
            {
                {1.0 * cH, -1.0 * cH },
                {-1.0 * cH, 1.0 * cH }
            };
            return stiffMatrix;
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double cH = CalculateConductivity();
            double gH = CalculateTemperatureJump();
            double[] heatFlux = new double[]
            {
                -1.0 * cH * gH ,
                 1.0 * cH * gH
            };
            return heatFlux;
        }

        public double[,] CreateMassMatrix()
        {
            return new double[4, 4];
        }

        public double[,] CreateDampingMatrix()
        {
            return new double[4, 4];
        }
    }
}

