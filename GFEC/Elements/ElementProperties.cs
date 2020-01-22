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

        public ElementProperties(double youngMod, double sectionArea, string elementType)
        {
            YoungMod = youngMod;
            SectionArea = sectionArea;
            ElementType = elementType;
        }

        public ElementProperties(double youngMod, double sectionArea, double momentOfInertia, string elementType)
        {
            YoungMod = youngMod;
            SectionArea = sectionArea;
            MomentOfInertia = momentOfInertia;
            ElementType = elementType;
        }

        public ElementProperties()
        {

        }
    }
}
