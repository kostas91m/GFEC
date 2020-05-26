using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public interface ISolver
    {
        ILinearSolution LinearScheme { get; set; }
        IAssembly AssemblyData { get; set; }
        void Solve(double[] rhsVector);
        bool ActivateNonLinearSolver { get; set; }
        INonLinearSolution NonLinearScheme { get; set; }
        void PrintSolution();
        double[,] CustomStiffnessMatrix { get; set; }
        double[] GetSolution();
        Dictionary<int, double[]> GetInternalForces();
        Dictionary<int, double[]> GetAllStepsSolutions();
    }
}