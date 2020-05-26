using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GFEC
{
    public static class AssemblyHelpMethods
    {
        public static Dictionary<int,double> RetrieveContactContactivity(IAssembly assembly)
        {
            Dictionary<int, double> contactElementContactivity = new Dictionary<int, double>();
            foreach (KeyValuePair<int, IElement> element in assembly.ElementsAssembly)
            {
                if (element.Value.GetType() == typeof(ContactNtN2DTh))
                {
                    double contactivity = ((ContactNtN2DTh)element.Value).CalculateConductivity();
                    int elementID = element.Key;
                    contactElementContactivity.Add(elementID, contactivity);
                }
            }
            return contactElementContactivity;
        }
    }
}
