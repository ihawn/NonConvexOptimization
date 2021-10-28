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
            for(int j = i + 1 ; j < hypArr.Length; j++)
            {
                for(int k = j + 1; k < hypArr.Length; k++)
                {
                    if (hypArr[i].direction != hypArr[j].direction &&
                        hypArr[j].direction != hypArr[k].direction &&
                        hypArr[i].direction != hypArr[k].direction)
                    {
                        Vector3 pt = IntersectHyperplanes(hypArr[i], hypArr[j], hypArr[k]);
                        
                        if(ValidIntersection(p1, p2, p3, pt))
                            intPoints.Add(pt);
                    }
                }
            }
        }

        return intPoints;
    }

    public List<Vector3> IntersectNew(List<Hyperplane> lst, List<Pyramid> pyrLst, List<Vector3> sectList, Pyramid pyr)
    {
        for(int i = 0; i < lst.Count; i++)
        {
            for(int j = i+1; j < lst.Count; j++)
            {
                if(lst[i].direction != lst[j].direction)
                {
                    for (int k = 0; k < pyr.hyperplanes.Length; k++)
                    {
                        if (lst[j].direction != pyr.hyperplanes[k].direction &&
                            lst[i].direction != pyr.hyperplanes[k].direction)
                        {
                            Vector3 pt = IntersectHyperplanes(lst[i], lst[j], pyr.hyperplanes[k]);

                            if (ValidIntersection(pyrLst[lst[i].parentID], pyrLst[lst[j].parentID], pyr, pt))
                                sectList.Add(pt);
                        }
                    }
                }
            }
        }

        return sectList;
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

    //Checks whether an intersection is on p1, p2, and p3 or not
    bool ValidIntersection(Pyramid p1, Pyramid p2, Pyramid p3, Vector3 pt)
    {
        // y - L||x - hatx||
        float t1 = p1.peak.y - p1.L * Mathf.Max(Mathf.Abs(pt.x - p1.peak.x), Mathf.Abs(pt.z - p1.peak.z));
        float t2 = p2.peak.y - p2.L * Mathf.Max(Mathf.Abs(pt.x - p2.peak.x), Mathf.Abs(pt.z - p2.peak.z));
        float t3 = p3.peak.y - p3.L * Mathf.Max(Mathf.Abs(pt.x - p3.peak.x), Mathf.Abs(pt.z - p3.peak.z));

        return Mathf.Abs(pt.y - t1) <= 0.0000001f &&
               Mathf.Abs(pt.y - t2) <= 0.0000001f &&
               Mathf.Abs(pt.y - t3) <= 0.0000001f &&
               Vector3.Magnitude(pt) < 100000;
    }


    //Removes invalid intersections based on all pyramids, not just those which generated an intersection
    public List<Vector3> PruneIntersections(List<Vector3> sects, List<Pyramid> pyrs)
    {
        for(int i = 0; i < pyrs.Count; i++)
        {
            int s = sects.Count;
            int j = 0;
            while(j < s)
            {
                float t = pyrs[i].peak.y - pyrs[i].L * Mathf.Max(Mathf.Abs(sects[j].x - pyrs[i].peak.x), Mathf.Abs(sects[j].z - pyrs[i].peak.z));
                
                if(t > sects[j].y + 0.000001f)
                {
                    sects.RemoveAt(j);
                    s--;
                }

                j++;
            } 
        }

        return sects;
    }
}
