using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public static class AssemblyHelpMethods
    {
        public static List<double> RetrieveContactContactivity(IAssembly assembly)
        {
            List<double> contactElementContactivity = new List<double>();
            foreach (KeyValuePair<int, IElement> element in assembly.ElementsAssembly)
            {
                if (element.Value.GetType() == typeof(ContactNtN2DTh))
                {
                    double contactivity = ((ContactNtN2DTh)element.Value).CalculateConductivity();
                    contactElementContactivity.Add(contactivity);
                }
            }
            return contactElementContactivity;
        }
    }
}
