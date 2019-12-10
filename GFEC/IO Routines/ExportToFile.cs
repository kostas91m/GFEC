using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace GFEC
{
    public static class ExportToFile
    {
        public static void ExportExplicitResults(Dictionary<int, double[]> solution, Dictionary<int, double> timeAtEachStep, int dofNumber, int intervals)
        {
            string[] lines = new string[solution.Count / intervals];
            int step = 0;
            int line = 0;
            while (step < solution.Count - 1)
            {
                double[] sol = solution[step];
                lines[line] = line.ToString() + " " + timeAtEachStep[step].ToString() + " " + sol[dofNumber].ToString();
                line = line + 1;
                step = step + intervals;
                //if (step >= solution.Count-1)
                //{
                //    break;
                //}
            }
            //File.WriteAllLines(@"D:\WriteLines2.txt", lines);
        }

        public static void ExportGeometryDataWithTemperatures(Dictionary<int, INode> nodesList, double[] temperatures)
        {
            if (nodesList.Count != temperatures.Length)
            {
                throw new Exception("Mismatch beetween total nodes and temperatures vector");
            }

            string[] lines = new string[nodesList.Count];
            for (int i = 1; i <= nodesList.Count; i++)
            {
                lines[i - 1] = nodesList[i].XCoordinate.ToString() + " " + nodesList[i].YCoordinate.ToString() + " " + temperatures[i - 1];
            }
            File.WriteAllLines(@"C:\Results.dat", lines);
        }
    }
}
