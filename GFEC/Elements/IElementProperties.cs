using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public interface IElementProperties
    {
        double YoungMod { get; set; }
        double SectionArea { get; set; }
        double MomentOfInertia { get; set; }
        string ElementType { get; set; }
        double Density { get; set; }
        double Thickness { get; set; }
        double ThermalConductivity { get; set; }
        double ContactForceValue { get; set; }
        double ContactThermalConductivity { get; set; }
        double SurfaceRoughness { get; set; }
        double YieldStrength { get; set; }
        double Dx1 { get; set; }
        double Dx2 { get; set; }
        double Dx { get; set; }
        double A { get; set; }
        double B { get; set; }
        double PoissonRatio { get; set; }

    }
}

