using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GameManager : MonoBehaviour
{
    public PyramidHandler ph;
    public MGrapher mg;
    public DotGrapher dg;
    public Intersections ih;
    
    void Start()
    {
        UnityEditor.SceneView.FocusWindowIfItsOpen(typeof(UnityEditor.SceneView));

        float L = 3;

        Pyramid pyr1 = ph.GeneratePyramid(new Vector3(0, 1, 0), L, 0);
        Pyramid pyr2 = ph.GeneratePyramid(new Vector3(0.3f, 1, 0.2f), L, 1);
        Pyramid pyr3 = ph.GeneratePyramid(new Vector3(-1.2f, 1.3f, 0.3f), L, 2);

        mg.GeneratePyramidMesh(pyr1);
        mg.GeneratePyramidMesh(pyr2);
        mg.GeneratePyramidMesh(pyr3);


        List<Vector3> intersections = ih.IntersectPyramids(pyr1, pyr2, pyr3);

        for(int i = 0; i < intersections.Count; i++)
        {
            dg.GraphPoint(intersections[i], 0.15f, Color.red);
        }
    }
}
