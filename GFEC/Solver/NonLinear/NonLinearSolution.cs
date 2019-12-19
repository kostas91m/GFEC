using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public class NonLinearSolution : INonLinearSolution
    {
        public int numberOfLoadSteps { get; set; } = 10;
        protected int[] boundaryDof;
        protected IAssembly discretization;
        protected double lambda;
        protected double tolerance = 1e-5;
        protected int maxIterations = 100;
        public bool PrintResidual { get; set; } = false;
        protected ILinearSolution linearSolver;
        public Dictionary<int, double[]> InternalForces { get; set; }
        public Dictionary<int, double[]> Solutions { get; set; }

        public virtual double[] Solve(IAssembly assembly, ILinearSolution linearScheme, double[] forceVector)
        {
            throw new Exception("LinearSolution.Solve not implemented");
        }


    }
}
