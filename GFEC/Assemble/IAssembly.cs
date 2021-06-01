using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GFEC
{
    public interface IAssembly
    {
        Dictionary<int, IElementProperties> ElementsProperties { get; set; }
        Dictionary<int, INode> Nodes { get; set; }
        Dictionary<int, Dictionary<int, int>> ElementsConnectivity { get; set; }
        Dictionary<int, bool[]> NodeFreedomAllocationList { get; set; }
        void CreateElementsAssembly();
        double[,] CreateTotalStiffnessMatrix();
        bool ActivateBoundaryConditions { get; set; }
        int[] BoundedDOFsVector { get; set; }
        void UpdateDisplacements(double[] totalDisplacementVector);
        void UpdateAccelerations(double[] totalAccelerationsVector);
        double[] CreateTotalInternalForcesVector();
        double[,] CreateTotalMassMatrix();
        double[,] CreateTotalDampingMatrix();
        Dictionary<int, double[]> GetElementsInternalForces(double[] totalInternalForcesVector);
        Dictionary<int, List<double[]>> GetElementsStresses(double[] totalDisplacementVector);
        Dictionary<int, List<double[]>> GetElementsStains(double[] totalDisplacementVector);
        Dictionary<int, List<double[]>> GetElementsNodesStresses(double[] totalDisplacementVector);
        Dictionary<int, List<double[]>> GetElementsNodesStains(double[] totalDisplacementVector);
        Dictionary<int, List<double[]>> GetElementsGaussPoints(double[] totalDisplacementVector);


        List<string> GetElementsType();
        Dictionary<int, IElement> ElementsAssembly { get; set; }
        int CountElementsOfSameType(Type elementType);




        //void UpdateValues(double[] totalDisplacementVector);
        //double[,] CreateTotalStiffnessMatrix();
        //double[,] CreateTotalMassMatrix();
        //double[] CreateTotalInternalForcesVector();
        //int[] BoundedDOFsVector
        //{
        //    get;
        //    set;
        //}
    }
}