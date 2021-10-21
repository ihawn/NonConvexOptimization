using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;

public class Intersections : MonoBehaviour
{
    public List<Vector3> IntersectPyramids(Pyramid p1, Pyramid p2, Pyramid p3)
    {
        Hyperplane[] hypArr = new Hyperplane[12];
        for(int i = 0; i < 4; i++)
        {
            hypArr[3*i] = p1.hyperplanes[i];
            hypArr[3*i + 1] = p2.hyperplanes[i];
            hypArr[3*i + 2] = p3.hyperplanes[i];
        }

        List<Vector3> intPoints = new List<Vector3>();

        for(int i = 0; i < hypArr.Length; i++)
        {
            for(int j = i +1 ; j < hypArr.Length; j++)
            {
                for(int k = j + 1; k < hypArr.Length; k++)
                {
                    var A = Matrix<double>.Build.DenseOfArray(new double[,]
                    {
                        { hypArr[i].coeff[0], hypArr[i].coeff[1], hypArr[i].coeff[2] },
                        { hypArr[j].coeff[0], hypArr[j].coeff[1], hypArr[j].coeff[2] },
                        { hypArr[k].coeff[0], hypArr[k].coeff[1], hypArr[k].coeff[2] },
                    });

                    if (A.Determinant() != 0)
                    {
                        Vector3 pt = IntersectHyperplanes(hypArr[i], hypArr[j], hypArr[k]);
                        
                        if(ValidIntersection(p1, p2, p3, pt))
                        {
                            intPoints.Add(pt);
                        }
                    }
                }
            }
        }

        return intPoints;
    }

    Vector3 IntersectHyperplanes(Hyperplane h1, Hyperplane h2, Hyperplane h3)
    {
        var A = Matrix<double>.Build.DenseOfArray(new double[,]
        {
            { h1.coeff[0], h1.coeff[1], h1.coeff[2] },
            { h2.coeff[0], h2.coeff[1], h2.coeff[2] },
            { h3.coeff[0], h3.coeff[1], h3.coeff[2] },
        });
        var b = Vector<double>.Build.Dense(new double[] { -h1.coeff[3], -h2.coeff[3], -h3.coeff[3] });
        var x = A.Solve(b);

        return new Vector3((float) x[0], (float) x[1], (float) x[2]);
    }

    bool ValidIntersection(Pyramid p1, Pyramid p2, Pyramid p3, Vector3 pt)
    {
        float t1 = p1.peak.y - p1.L * Mathf.Max(Mathf.Abs(pt.x - p1.peak.x), Mathf.Abs(pt.z - p1.peak.z));
        float t2 = p2.peak.y - p2.L * Mathf.Max(Mathf.Abs(pt.x - p2.peak.x), Mathf.Abs(pt.z - p2.peak.z));
        float t3 = p3.peak.y - p3.L * Mathf.Max(Mathf.Abs(pt.x - p3.peak.x), Mathf.Abs(pt.z - p3.peak.z));

        return Mathf.Abs(pt.y - t1) <= 0.00001f &&
               Mathf.Abs(pt.y - t2) <= 0.00001f &&
               Mathf.Abs(pt.y - t3) <= 0.00001f &&
               Vector3.Magnitude(pt) < 100000;
    }
}
