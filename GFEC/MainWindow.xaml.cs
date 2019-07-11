using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using LiveCharts;
using LiveCharts.Wpf;

namespace GFEC
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public SeriesCollection Graph {get; set;}
        public SeriesCollection Mesh { get; set; }
        private Results solverResults;

        public MainWindow()
        {
            InitializeComponent();
            LoadComboBox();
            
        }

        private void RunButton(object sender, RoutedEventArgs args)
        {
            SolveSelectedExample();
            Graph = ShowToGUI.ShowResults(solverResults, 1, 1);
            

            Dictionary<int, INode> nodes = new Dictionary<int, INode>();
            nodes[1] = new Node(0.0, 0.01);
            nodes[2] = new Node(0.3, 0.01);
            nodes[3] = new Node(0.6, 0.01);
            nodes[4] = new Node(0.6, 0.12);
            nodes[5] = new Node(0.3, 0.12);
            nodes[6] = new Node(0.0, 0.12);
            nodes[7] = new Node(0.45, -0.11);
            nodes[8] = new Node(0.75, -0.11);
            nodes[9] = new Node(1.05, -0.11);
            nodes[10] = new Node(0.45, 0.0);
            nodes[11] = new Node(0.75, 0.0);
            nodes[12] = new Node(1.05, 0.0);
            Dictionary<int, Dictionary<int, int>> connectivity = new Dictionary<int, Dictionary<int, int>>();
            connectivity[1] = new Dictionary<int, int>() { { 1, 1 }, { 2, 2 }, { 3, 5 }, { 4, 6 } };
            connectivity[2] = new Dictionary<int, int>() { { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 5 } };
            connectivity[3] = new Dictionary<int, int>() { { 1, 7 }, { 2, 8 }, { 3, 11 }, { 4, 10 } };
            connectivity[4] = new Dictionary<int, int>() { { 1, 8 }, { 2, 9 }, { 3, 12 }, { 4, 11 } };
            //connectivity[5] = new Dictionary<int, int>() { { 1, 10 }, { 2, 11 }, { 3, 3 } };
            Mesh = ShowToGUI.DrawMesh(nodes, connectivity);

            DataContext = this;

            return;
        }

        private void LoadComboBox()
        {
            List<string> exampleList = new List<string>();
            exampleList.Add("LinearTrussExample");
            exampleList.Add("TwoQuadsExample");
            exampleList.Add("TwoBeamsInFrContactQuadsExample");

            ComboBox1.ItemsSource = exampleList;
        }

        private void SolveSelectedExample()
        {
            Results finalResults;
            Tuple<Dictionary<int, double[]>, Dictionary<int, double>> results;
              string selectedExample = ComboBox1.SelectedItem.ToString();
            switch (selectedExample)
            {
                case "TwoQuadsExample":
                    finalResults = TwoQuadsExample.RunStaticExample();
                    break;
                case "LinearTrussExample":
                    finalResults = LinearTrussExample.RunExample();
                    break;
                case "TwoBeamsInFrContactQuadsExample":
                    finalResults = TwoBeamsInFrContactQuadsExample.RunDynamicExample();
                    break;
                default:
                    finalResults = TwoQuadsExample.RunStaticExample();
                    break;
            }
            //Results.Text = solution[0].ToString();
            
            solverResults = finalResults;
        }
    }

    
}
