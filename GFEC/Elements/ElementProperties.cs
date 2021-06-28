using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public class ElementProperties : IElementProperties
    {
        public double YoungMod { get; set; }
        public double SectionArea { get; set; }
        public double MomentOfInertia { get; set; }
        public string ElementType { get; set; }
        public double Density { get; set; }
        public double Thickness { get; set; }
        public double ThermalConductivity { get; set; }
        public double ContactForceValue { get; set; }

        public double ContactThermalConductivity { get; set; }
        public double SurfaceRoughness { get; set; }
        public double YieldStrength { get; set; }
        public double Dx1 { get; set; }
        public double Dx2 { get; set; }
        public double Dx { get; set; }
        public double A { get; set; }
        public double B { get; set; }
        public double PoissonRatio { get; set; }
        public int MasterSegmentPolynomialDegree { get; set; }
        public int SlaveSegmentPolynomialDegree { get; set; }
        public int IntegrationPoints { get; set; }
        public double PenaltyFactorRatio { get; set; }

        public ElementProperties(double youngMod, double sectionArea, string elementType)
        {
            YoungMod = youngMod;
            SectionArea = sectionArea;
            ElementType = elementType;
            PoissonRatio = 0.30;
            PenaltyFactorRatio = 10.0;
        }
        /// <summary>
        /// New constructor after the inclusion of Poissons ratio in element properties
        /// </summary>
        /// <param name="youngMod"></param>
        /// <param name="poissonRatio"></param>
        /// <param name="sectionArea"></param>
        /// <param name="thickness"></param>
        /// <param name="density"></param>
        /// <param name="elementType"></param>
        public ElementProperties(double youngMod, double poissonRatio, double sectionArea, double thickness, double density, string elementType)
        {
            YoungMod = youngMod;
            SectionArea = sectionArea;
            ElementType = elementType;
            PoissonRatio = poissonRatio;
            Thickness = thickness;
            Density = density;
            PenaltyFactorRatio = 10.0;
        }
        public ElementProperties(double youngMod, double sectionArea, double momentOfInertia, string elementType)
        {
            YoungMod = youngMod;
            SectionArea = sectionArea;
            MomentOfInertia = momentOfInertia;
            ElementType = elementType;
            PoissonRatio = 0.30;
            PenaltyFactorRatio = 10.0;
        }
        /// <summary>
        /// New constructor for STS contact element. Properties are added in order to 
        /// include higher order elements
        /// </summary>
        public ElementProperties(double youngMod, double contactArea, string elementType, double penaltyFactorRatio,
            int integrationPoints, int slaveSegmentPolynomialDegree, int masterSegmentPolynomialDegree)
        {
            if (youngMod > 0.0)
            {
                YoungMod = youngMod;

            }
            else
            {
                throw new Exception("Young's modulus must be positive real number");

            }
            SectionArea = contactArea;
            if (penaltyFactorRatio > 0.0)
            {
                PenaltyFactorRatio = penaltyFactorRatio;

            }
            else
            {
                PenaltyFactorRatio = 10.0;//default value

            }
            ElementType = elementType;
            if (integrationPoints <= 10 && integrationPoints >= 1)
            {
                IntegrationPoints = integrationPoints;

            }
            else
            {
                throw new Exception("The amount of integration points must be between 1 & 10");
            }
            if (slaveSegmentPolynomialDegree<=2 && slaveSegmentPolynomialDegree >= 1)
            {
                SlaveSegmentPolynomialDegree = slaveSegmentPolynomialDegree;

            }
            else
            {
                throw new Exception("Only polynomial degrees <= 2");
            }
            if (masterSegmentPolynomialDegree <= 2 && masterSegmentPolynomialDegree >= 1)
            {
                MasterSegmentPolynomialDegree = masterSegmentPolynomialDegree;

            }
            else
            {
                throw new Exception("Only polynomial degrees <= 2");
            }
        }
        public ElementProperties()
        {

        }
    }
}
