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
    public class ShowDiagramInGUIArgs: EventArgs
    {
        public SeriesCollection DiagramData { get; set; }
    }
}
