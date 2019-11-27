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