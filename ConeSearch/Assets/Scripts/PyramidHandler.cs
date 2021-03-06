using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;

public class PyramidHandler : MonoBehaviour
{
    public GameManager gm;
    public DotGrapher dg;

    public Hyperplane GenerateHyperplane(float[,] x, int parID, int direction)
    {
        var M = Matrix<float>.Build;
        var V = Vector<float>.Build;
        var m = M.DenseOfArray(new[,] { { x[0,0],  x[0,1], x[0,2], 1.0f },
                                        { x[1,0],  x[1,1], x[1,2], 1.0f },
                                        { x[2,0],  x[2,1], x[2,2], 1.0f },
                                        { 1.0f, 1.0f, 1.0f, 1.0f} });
        
        var b = Vector<float>.Build.Dense(new float[] { 0f, 0f, 0f, 1f });
        var z = m.Solve(b);
        var c = V.DenseOfArray(new[] { z[0], z[1], z[2], z[3] });

        return new Hyperplane(parID, direction, c);
    }

    public Pyramid GeneratePyramid(Vector3 peak, float L, int id)
    {
        peak.y = gm.Objective(peak.x, peak.z);

        float[,] x1 = new float[,]
        {
            { peak.x + 1, peak.y - L, peak.z + 1 },
            { peak.x - 1, peak.y - L, peak.z + 1 },
            { peak.x, peak.y, peak.z }
        };
        float[,] x2 = new float[,]
        {
            { peak.x - 1, peak.y - L, peak.z + 1 },
            { peak.x - 1, peak.y - L, peak.z - 1 },
            { peak.x, peak.y, peak.z }
        };
        float[,] x3 = new float[,]
        {
            { peak.x - 1, peak.y - L, peak.z - 1 },
            { peak.x + 1, peak.y - L, peak.z - 1 },
            { peak.x, peak.y, peak.z }
        };
        float[,] x4 = new float[,]
        {
            { peak.x + 1, peak.y - L, peak.z - 1 },
            { peak.x + 1, peak.y - L, peak.z + 1 },
            { peak.x, peak.y, peak.z }
        };

        Hyperplane[] hyps = new Hyperplane[]
        {
            GenerateHyperplane(x1, id, 0),
            GenerateHyperplane(x2, id, 1),
            GenerateHyperplane(x3, id, 2),
            GenerateHyperplane(x4, id, 3)
        };

        return new Pyramid(id, peak, L, 0, hyps);
    }

    public List<Pyramid[]> CombinePyramids(List<Pyramid> pyrList)
    {
        List<Pyramid[]> combos = new List<Pyramid[]>();

        for (int i = 0; i < pyrList.Count; i++)
        {
            for(int j = i+1; j < pyrList.Count; j++)
            {
                for(int k = j+1; k < pyrList.Count; k++)
                {
                    combos.Add(new Pyramid[] { pyrList[i], pyrList[j], pyrList[k] });
                }
            }
        }
        return combos;
    }
}

