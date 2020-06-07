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
                Type thermalContactElementType = element.Value.GetType();
                int elementID = element.Key;
                if (thermalContactElementType == typeof(ContactNtN2DTh))
                {
                    double contactivity = ((ContactNtN2DTh)element.Value).CalculateConductivity();
                    //int elementID = element.Key;
                    contactElementContactivity.Add(elementID, contactivity);
                }
                else if (thermalContactElementType == typeof(ContactNtS2DTh))
                {
                    double contactivity = ((ContactNtS2DTh)element.Value).CalculateConductivity();
                    //int elementID = element.Key;
                    contactElementContactivity.Add(elementID, contactivity);
                }
            }
            return contactElementContactivity;
        }
    }
}
