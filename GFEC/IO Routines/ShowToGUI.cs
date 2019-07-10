using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using LiveCharts;
using LiveCharts.Configurations;
using LiveCharts.Defaults;
using LiveCharts.Wpf;

namespace GFEC
{
    public static class ShowToGUI
    {
        public static SeriesCollection ShowResults(Results analysisResults, int dofNumber, int intervals)
        {
            int countVector = analysisResults.DynamicSolution.Count;
            int step = 0;
            int line = 0;
            double[] xAxis = new double[analysisResults.DynamicSolution.Count];
            double[] yAxis = new double[analysisResults.DynamicSolution.Count];
            var points = new ObservablePoint[analysisResults.DynamicSolution.Count];
            while (step < analysisResults.DynamicSolution.Count)
            {
                double[] sol = analysisResults.DynamicSolution[step];
                xAxis[line] = analysisResults.TimeSteps[step];
                yAxis[line] = sol[dofNumber];
                points[line] = new ObservablePoint() { X = analysisResults.TimeSteps[step], Y = sol[dofNumber] };
                line = line + 1;
                step = step + intervals;
                
                //if (step >= solution.Count-1)
                //{
                //    break;
                //}
            }
            var mapper = Mappers.Xy<ObservablePoint>() //in this case value is of type <ObservablePoint>
                .X(value => value.X) //use the X property as X
                .Y(value => value.Y); //use the Y property as Y
            SeriesCollection graph = new SeriesCollection
            {
                //new LineSeries
                //{
                //    Values = new ChartValues<double>(xAxis)
                //},
                new LineSeries
                {
                    Values = new ChartValues<ObservablePoint>(points)
                }
            };
            return graph;
        }
    }
}
