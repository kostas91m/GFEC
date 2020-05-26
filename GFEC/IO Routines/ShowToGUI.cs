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
using System.Diagnostics;
using System.IO;
using System.Threading;
using System.Windows;

namespace GFEC
{
    public class ShowToGUI
    {
        public event EventHandler<ShowDiagramInGUIArgs> ShowDiagramInGUI;
        public event EventHandler<SeriesCollection> TestEvent;

        protected virtual void OnShowDiagramInGUI(ShowDiagramInGUIArgs e)
        {
            EventHandler<ShowDiagramInGUIArgs> handler = ShowDiagramInGUI;
            if (handler != null)
            {
                handler(this, e);
            }
        }

        protected virtual void OnTestEvent(SeriesCollection e)
        {
            EventHandler<SeriesCollection> handler = TestEvent;
            if (handler!=null)
            {
                handler(this, e);
            }
        }

        public void TestEventMethod()
        {
            SeriesCollection Something = new SeriesCollection
            {
                new LineSeries
                {
                    Values = new ChartValues<double> { 3, 5, 7, 4 }
                },
                new ColumnSeries
                {
                    Values = new ChartValues<decimal> { 5, 6, 2, 7 }
                }
            };

            OnTestEvent(Something);
        }

        public SeriesCollection ShowResults(Results analysisResults)
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
                    return ShowStaticNonLinearResults(analysisResults);
                    break;
            }
        }

        private SeriesCollection ShowDynamicLinearResults(Results analysisResults)
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

        private SeriesCollection ShowStaticNonLinearResults(Results analysisResults)
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

            OnShowDiagramInGUI(new ShowDiagramInGUIArgs(){ DiagramData = graph});
            return graph;
        }

        public static void PlotInitialGeometry(IAssembly assembly)
        {
            GnuPlot.Set("terminal png size 1920, 1080");
            GnuPlot.Set("output 'gnuplot.png'");
            GnuPlot.HoldOn();
            GnuPlot.Set("size ratio -1");
            GnuPlot.Unset("key");
            double[] X,Y;
            List<StoredPlot> storedPlots = new List<StoredPlot>();
            //List<double> X = new List<double>();
            //List<double> Y = new List<double>();
            foreach (var element in assembly.ElementsAssembly)
            {    
                int numberOfNodes = element.Value.Nodes.Count;
                for (int i = 1; i <= numberOfNodes; i++)
                {
                    if (i != numberOfNodes)
                    {
                        //X.Add(element.Value.Nodes[i].XCoordinate);
                        //X.Add(element.Value.Nodes[i+1].XCoordinate);
                        //Y.Add(element.Value.Nodes[i].YCoordinate);
                        //Y.Add(element.Value.Nodes[i + 1].YCoordinate);
                        X = new double[] { element.Value.Nodes[i].XCoordinate, element.Value.Nodes[i + 1].XCoordinate };
                        Y = new double[] { element.Value.Nodes[i].YCoordinate, element.Value.Nodes[i + 1].YCoordinate };
                    }
                    else
                    {
                        //X.Add(element.Value.Nodes[i].XCoordinate);
                        //X.Add(element.Value.Nodes[ 1].XCoordinate);
                        //Y.Add(element.Value.Nodes[i].YCoordinate);
                        //Y.Add(element.Value.Nodes[ 1].YCoordinate);
                        X = new double[] { element.Value.Nodes[i].XCoordinate, element.Value.Nodes[1].XCoordinate };
                        Y = new double[] { element.Value.Nodes[i].YCoordinate, element.Value.Nodes[1].YCoordinate };
                    }
                    
                    if (element.Value is ContactNtN2D)
                    {
                        storedPlots.Add(new StoredPlot(X, Y, "with linespoints pt " + (int)PointStyles.SolidCircle + " lt rgb \"red\""));
                        //GnuPlot.Plot(X, Y, "with linespoints pt " + (int)PointStyles.SolidCircle + " lt rgb \"blue\"");
                    }
                    else if (element.Value is ContactNtS2D)
                    {
                        storedPlots.Add(new StoredPlot(X, Y, "with linespoints pt " + (int)PointStyles.SolidCircle + " lt rgb \"red\""));
                    }
                    else
                    {
                        storedPlots.Add(new StoredPlot(X, Y, "with linespoints pt " + (int)PointStyles.SolidCircle + " lt rgb \"blue\""));
                        //GnuPlot.Plot(X, Y, "with linespoints pt " + (int)PointStyles.SolidCircle + " lt rgb \"blue\"");
                    }


                }
            }
            GnuPlot.Plot(storedPlots);
            GnuPlot.HoldOff();
            //GnuPlot.Close();
            //GnuPlot.KillProcess();

        }

        public static void PlotFinalGeometry(IAssembly assembly)
        {
            GnuPlot.Set("terminal png size 1920, 1080");
            GnuPlot.Set("output 'gnuplot2.png'");
            GnuPlot.HoldOn();
            GnuPlot.Set("size ratio -1");
            GnuPlot.Unset("key");
            double[] X, Y;
            List<StoredPlot> storedPlots = new List<StoredPlot>();
            
            foreach (var element in assembly.ElementsAssembly)
            {
                //int numberOfNodes = element.Value.Nodes.Count;
                Dictionary<int, INode> nodesAtFinalState = element.Value.NodesAtFinalState();
                int numberOfNodes = nodesAtFinalState.Count;
                for (int i = 1; i <= numberOfNodes; i++)
                {
                    if (i != numberOfNodes)
                    {
                        X = new double[] { nodesAtFinalState[i].XCoordinate, nodesAtFinalState[i + 1].XCoordinate };
                        Y = new double[] { nodesAtFinalState[i].YCoordinate, nodesAtFinalState[i + 1].YCoordinate };
                    }
                    else
                    {
                        X = new double[] { nodesAtFinalState[i].XCoordinate, nodesAtFinalState[1].XCoordinate };
                        Y = new double[] { nodesAtFinalState[i].YCoordinate, nodesAtFinalState[1].YCoordinate };
                    }

                    if (element.Value is ContactNtN2D)
                    {
                        storedPlots.Add(new StoredPlot(X, Y, "with linespoints pt " + (int)PointStyles.SolidCircle + " lt rgb \"red\""));
                    }
                    else
                    {
                        storedPlots.Add(new StoredPlot(X, Y, "with linespoints pt " + (int)PointStyles.SolidCircle + " lt rgb \"blue\""));
                    }


                }
            }
            GnuPlot.Plot(storedPlots);
            GnuPlot.HoldOff();
            //GnuPlot.Close();
            //GnuPlot.KillProcess();

        }

        public static void PlotHeatMap(List<HeatMapData> plots)
        {
            GnuPlot.HoldOn();
            GnuPlot.Set("cbrange[0:7.0]");
            GnuPlot.Set("palette defined(0 \"blue\", 0.33\"green\", 0.67\"yellow\", 1 \"red\")");
            GnuPlot.Set("pm3d");
            GnuPlot.Set("dgrid3d");
            GnuPlot.Set("view map");
            foreach (var plot in plots)
            {
                GnuPlot.SPlot(plot.Xcoordinates, plot.Ycoordinates, plot.Temperatures);
            }

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
