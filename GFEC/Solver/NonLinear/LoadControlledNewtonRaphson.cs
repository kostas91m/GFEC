using System;
using System.Collections.Generic;
using System.Threading;

namespace GFEC
{
    public class LoadControlledNewtonRaphson : NonLinearSolution
    {
        private double[] localSolutionVector;

        public LoadControlledNewtonRaphson()
        {

        }
        public LoadControlledNewtonRaphson(double[] exSolution)
        {
            localSolutionVector = exSolution;
        }
        private double[] LoadControlledNR(double[] forceVector)
        {
            double[] incrementDf = VectorOperations.VectorScalarProductNew(forceVector, lambda);
            double[] solutionVector = localSolutionVector;
            double[] incrementalExternalForcesVector = new double[forceVector.Length];
            double[] tempSolutionVector = new double[solutionVector.Length];
            double[] deltaU = new double[solutionVector.Length];
            double[] internalForcesTotalVector;
            double[] dU;
            double[] residual;
            double residualNorm;
            for (int i = 0; i < numberOfLoadSteps; i++)
            {
                incrementalExternalForcesVector = VectorOperations.VectorVectorAddition(incrementalExternalForcesVector, incrementDf);
                discretization.UpdateDisplacements(solutionVector);
                internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();
                double[,] stiffnessMatrix = discretization.CreateTotalStiffnessMatrix();
                //OnConvergenceResult("Newton-Raphson: Solution not converged at load step" + i); 
                dU = linearSolver.Solve(stiffnessMatrix, incrementDf);
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, dU);
                residual = VectorOperations.VectorVectorSubtraction(internalForcesTotalVector, incrementalExternalForcesVector);
                residualNorm = VectorOperations.VectorNorm2(residual);
                int iteration = 0;
                Array.Clear(deltaU, 0, deltaU.Length);
                while (residualNorm > Tolerance && iteration < MaxIterations)
                {
                    stiffnessMatrix = discretization.CreateTotalStiffnessMatrix();
                    deltaU = VectorOperations.VectorVectorSubtraction(deltaU, linearSolver.Solve(stiffnessMatrix, residual));
                    tempSolutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                    discretization.UpdateDisplacements(tempSolutionVector);
                    internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();
                    residual = VectorOperations.VectorVectorSubtraction(internalForcesTotalVector, incrementalExternalForcesVector);
                    residualNorm = VectorOperations.VectorNorm2(residual);
                    if (residualNorm <= Tolerance)
                    {
                        OnConvergenceResult("Newton-Raphson: Load Step "+i+ " - Solution converged at iteration " + iteration + " - Residual Norm = " +residualNorm);
                    }
                    else
                    {
                        OnConvergenceResult("Newton-Raphson: Load Step " + i + " - Solution not converged at iteration " + iteration + " - Residual Norm = " + residualNorm);
                    }
                    iteration = iteration + 1;
                    //(Application.Current.Windows[0] as MainWindow).LogTool.Text = "ok"; 
                    //OnConvergenceResult("Newton-Raphson: Solution not converged at load step" + iteration);
                }
                InternalForces.Add(i + 1, internalForcesTotalVector);
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                Solutions.Add(i + 1, solutionVector);
                if (iteration >= MaxIterations)
                {
                    OnConvergenceResult("Newton-Raphson did not converge at Load Step " + i + ". Exiting solution.");
                    LoadStepConvergence.Add("Solution not converged.");
                    break;
                }
                LoadStepConvergence.Add("Solution converged.");

            }
            return solutionVector;
        }
        
        public override double[] Solve(IAssembly assembly, ILinearSolution linearScheme, double[] forceVector)
        {
            InternalForces = new Dictionary<int, double[]>();
            Solutions = new Dictionary<int, double[]>();
            LoadStepConvergence = new List<string>();
            if (localSolutionVector == null)
            {
                localSolutionVector = new double[forceVector.Length];
            }
            discretization = assembly;
            linearSolver = linearScheme;
            lambda = 1.0 / numberOfLoadSteps;
            //double[] solution = null;

            //Thread tcore1 = new Thread(() => LoadControlledNR(forceVector));
            //tcore1.Start();
            //tcore1.Join();
            double[] solution = LoadControlledNR(forceVector);
            return solution;
        }

    }
}

