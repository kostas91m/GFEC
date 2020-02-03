using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    interface INonLinearSolution
    {
        double[] Solve(IAssembly assembly, ILinearSolution linearScheme, double[] forceVector);
        int numberOfLoadSteps { get; set; }
        Dictionary<int, double[]> InternalForces { get; set; }
        Dictionary<int, double[]> Solutions { get; set; }
        double Tolerance { get; set; }
        int MaxIterations { get; set; }
    }
}