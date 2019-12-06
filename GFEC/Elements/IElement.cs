using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public interface IElement
    {
        Dictionary<int, INode> Nodes { get; }
        IElementProperties Properties { get; set; }
        Dictionary<int, bool[]> ElementFreedomSignature { get; }
        List<int> ElementFreedomList { get; set; }
        //double CalculateElementLength();
        //double CalculateElementSinus();
        //double CalculateElementCosinus();
        //double[,] CreateLocalStiffnessMatrix();
        double[,] CreateGlobalStiffnessMatrix();
        double[,] CreateMassMatrix();
        double[] DisplacementVector { get; set; }
        double[] AccelerationVector { get; set; }
        double[] CreateInternalGlobalForcesVector();
        double[,] CreateDampingMatrix();
        Dictionary<int, INode> NodesAtFinalState();
    }
}

