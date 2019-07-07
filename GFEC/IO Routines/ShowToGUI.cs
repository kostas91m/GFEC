using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using LiveCharts;
using LiveCharts.Wpf;

namespace GFEC
{
    public static class ShowToGUI
    {
        public static SeriesCollection ShowResults(Results analysisResults, int dofNumber, int intervals)
        {
            int step = 0;
            int line = 0;
            double[] xAxis = new double[analysisResults.DynamicSolution.Count];
            double[] yAxis = new double[analysisResults.DynamicSolution.Count];
            while (step < analysisResults.DynamicSolution.Count - 1)
            {
                double[] sol = analysisResults.DynamicSolution[step];
                xAxis[line] = analysisResults.TimeSteps[step];
                yAxis[line] = sol[dofNumber];
                line = line + 1;
                step = step + intervals;
                //if (step >= solution.Count-1)
                //{
                //    break;
                //}
            }
            SeriesCollection graph = new SeriesCollection
            {
                new LineSeries
                {
                    Values = new ChartValues<double>(xAxis)
                },
                new LineSeries
                {
                    Values = new ChartValues<double>(yAxis)
                }
            };
            return graph;
        }
    }
}
