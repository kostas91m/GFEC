using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public static class MultiThreadingExample
    {
        public static void RunStaticExample()
        {
            double[,] matrix1 = new double[2000,2000];
            double[,] matrix2 = new double[2000,2000];
            double[,] result1, result2;

            Stopwatch watch1 = Stopwatch.StartNew();
            result1 = MatrixOperations.MatrixAddition(matrix1, matrix2);
            long first = watch1.ElapsedMilliseconds;

            Stopwatch watch2 = Stopwatch.StartNew();
            MatrixOperations.MatrixAdditionParallel(matrix1, matrix2);
            long second = watch2.ElapsedMilliseconds;
            result2 = MatrixOperations.TempVariable;

            Stopwatch watch3 = Stopwatch.StartNew();
            MatrixOperations.MatrixAdditionParallel2(matrix1, matrix2);
            long third = watch3.ElapsedMilliseconds;
        }
    }
}
