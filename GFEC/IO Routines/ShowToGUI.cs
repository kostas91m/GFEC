using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;
using LiveCharts;
using LiveCharts.Configurations;
using LiveCharts.Defaults;
using LiveCharts.Wpf;

namespace GFEC
{
    public static class ShowToGUI
    {
        public static SeriesCollection ShowResults(Results analysisResults)
        {
            switch (analysisResults.SolutionType)
            {
                case "Dynamic":
                    return ShowDynamicLinearResults(analysisResults);
                    break;
                case "Nonlinear":
                    return ShowStaticNonLinearResults(analysisResults);
                    break;
                default:
                    return ShowDynamicLinearResults(analysisResults);
                    break;
            }
        }

        private static SeriesCollection ShowDynamicLinearResults(Results analysisResults)
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
                yAxis[line] = sol[analysisResults.SelectedDOF];
                points[line] = new ObservablePoint() { X = analysisResults.TimeSteps[step], Y = sol[analysisResults.SelectedDOF] };
                line = line + 1;
                step = step + analysisResults.SelectedInterval;
                
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

        private static SeriesCollection ShowStaticNonLinearResults(Results analysisResults)
        {
            var points = new ObservablePoint[analysisResults.NonlinearSolution.Count];

            for (int i = 0; i < analysisResults.NonlinearSolution.Count; i++)
            {
                points[i] = new ObservablePoint() { X = i, Y = analysisResults.NonlinearSolution[i][analysisResults.SelectedDOF] };
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

        //public static SeriesCollection DrawMesh(Dictionary<int, INode> nodes, Dictionary<int, Dictionary<int, int>> connectivity)
        //{
        //    SeriesCollection mesh = new SeriesCollection();
        //    for (int j = 1; j <= connectivity.Count; j++)
        //    {
        //        var points = new ObservablePoint[connectivity[j].Count + 1];
        //        for (int i = 0; i < connectivity[j].Count; i++)
        //        {
        //            int nodeIndex = connectivity[j][i + 1];
        //            double x = nodes[nodeIndex].XCoordinate;
        //            double y = nodes[nodeIndex].YCoordinate;
        //            points[i] = new ObservablePoint() { X = x, Y = y };
        //        }
        //        points[points.Count() - 1] = new ObservablePoint();
        //        points[points.Count() - 1] = points[0];
        //        var mapper = Mappers.Xy<ObservablePoint>() //in this case value is of type <ObservablePoint>
        //            .X(value => value.X) //use the X property as X
        //            .Y(value => value.Y); //use the Y property as Y

        //        mesh.Add(new StepLineSeries
        //        {
        //            Values = new ChartValues<ObservablePoint>(points)
        //        });
        //    }

        //{
        //    new StepLineSeries
        //    {
        //        Values = new ChartValues<ObservablePoint>(points)
        //    }
        //};
        //return mesh;
        //}
    }
}
