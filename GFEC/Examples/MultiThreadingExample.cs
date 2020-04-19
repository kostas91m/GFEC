using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace GFEC
{
    public class MultiThreadingExample
    {
        public event EventHandler<string> timeElapsed;

        protected void OnTimeElapsed(string time)
        {
            if (timeElapsed != null)
            {
                timeElapsed.Invoke(this, time);
            }
        }
        public void RunExample()
        {
            double[,] matrix1 = MatrixOperations.CreateRandomMatrix(2000, 2000);
            double[,] matrix2 = MatrixOperations.CreateRandomMatrix(2000, 2000);
            double[] vector1 = VectorOperations.CreateRandomVector(2000);
            double[,] result1, result2, result3;
            double[] result1b, result2b, result3b;
            double result1c, result2c, result3c;
            MatrixOperations.ParallelCalculations = false;
            Stopwatch watch1 = Stopwatch.StartNew();
            result1 = MatrixOperations.MatrixAddition(matrix1, matrix2);
            result1b = VectorOperations.MatrixVectorProduct(result1, vector1);
            result1c = VectorOperations.VectorNorm2(result1b);
            long first = watch1.ElapsedMilliseconds;

            MatrixOperations.ParallelCalculations = true;
            Stopwatch watch2= Stopwatch.StartNew();
            result2 = MatrixOperations.MatrixAddition(matrix1, matrix2);
            //result2 = MatrixOperations.TempVariable;
            result2b = VectorOperations.MatrixVectorProduct(result2, vector1);
            result2c = VectorOperations.VectorNorm2(result2b);
            long second = watch2.ElapsedMilliseconds;
            

            string timeForCalculations = "Elapsed time for single threaded operation: " + first.ToString() + " -Result is:" + result1c + "\n" + "Elapsed time for multithreaded operation (Parallel for): " + second.ToString() + " -Result is:" + result2c;
            OnTimeElapsed(timeForCalculations);
        }
    }
}
