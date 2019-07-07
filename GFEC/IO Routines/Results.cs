using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public struct Results
    {
        public Dictionary<int, double[]> DynamicSolution { get; set; }
        public Dictionary<int, double> TimeSteps { get; set; }
    }
}
