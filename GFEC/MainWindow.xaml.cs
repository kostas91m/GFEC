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
