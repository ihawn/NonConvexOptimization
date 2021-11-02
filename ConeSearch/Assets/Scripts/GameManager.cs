using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;
using System;

public class GameManager : MonoBehaviour
{
    public PyramidHandler ph;
    public MGrapher mg;
    public DotGrapher dg;
    public Intersections ih;
    public List<Pyramid[]> combos;
    public List<Vector3> intersections;
    public List<Hyperplane> hyperplanes;
    public List<Pyramid> pyramids;
    

    public GameObject pointContainer;
    public float L = 3;
    public float bounds = 4f;
    public int closenessThreshold = 8;
    public int maxIterations = 5;
    public int pyramidCount = 0;
    public int hyperplaneCount = 0;


    public float lowerBound = -100000, upperBound = 100000;

    void Start()
    {
        UnityEditor.SceneView.FocusWindowIfItsOpen(typeof(UnityEditor.SceneView));

        pointContainer = new GameObject();
        hyperplanes = new List<Hyperplane>();


        Pyramid pyr1 = ph.GeneratePyramid(new Vector3(-bounds, 0, -bounds), L, 0);
        Pyramid pyr2 = ph.GeneratePyramid(new Vector3(-bounds, 0, bounds), L, 1);
        Pyramid pyr3 = ph.GeneratePyramid(new Vector3(bounds, 0, -bounds), L, 2);
        Pyramid pyr4 = ph.GeneratePyramid(new Vector3(bounds, 0, bounds), L, 3);

        pyramids = new List<Pyramid> { pyr1, pyr2, pyr3, pyr4 }; //List init

        //Draw the mesh
        mg.GeneratePyramidMesh(pyr1);
        mg.GeneratePyramidMesh(pyr2);
        mg.GeneratePyramidMesh(pyr3);
        mg.GeneratePyramidMesh(pyr4);


        //Combinations of pyramids intersecting
        combos = ph.CombinePyramids(pyramids);
        intersections = new List<Vector3>();
        for (int i = 0; i < combos.Count; i++)
        {
            intersections.AddRange(ih.IntersectPyramids(combos[i][0], combos[i][1], combos[i][2])); //Combine the lists
        }

        //Add hyperplanes to list
        for (int i = 0; i < pyramids.Count; i++)
        {
            for(int j = 0; j < pyramids[i].hyperplanes.Length; j++)
            {
                hyperplanes.Add(pyramids[i].hyperplanes[j]);

                intersections.AddRange(ih.IntersectNew(hyperplanes, pyramids, pyramids[i]));
                intersections = intersections.Distinct().ToList();
                ih.PruneIntersections(intersections, pyramids); //Remove intersection points that lie below any of the pyramids
            }    
        }


        //Draw the intersections
        for (int i = 0; i < intersections.Count; i++)
        {
            dg.GraphPoint(pointContainer, intersections[i], 0.03f, Color.red);
        }

        dg.GraphPoint(pointContainer, Vector3.zero, 0.015f, Color.blue); //Solution
    }

    void Update()
    {
        if(pyramids.Count < maxIterations)
        {
            int minLoc = MinSect(intersections); //Grab the min position in list
            Vector3 s = intersections[minLoc];
            float fx = Objective(s.x, s.z);
            Pyramid pyr = ph.GeneratePyramid(new Vector3(s.x, fx, s.z), L, pyramids.Count); //Generate new pyramid
            mg.GeneratePyramidMesh(pyr); //Draw the pyramid mesh


            if (s.y > lowerBound)
                lowerBound = s.y;
            if (fx < upperBound)
                upperBound = fx;


            List<Hyperplane> closeHyps = GetAdjHyperplanes(pyramids, pyr, closenessThreshold);
            List<Vector3> newIntersections = ih.IntersectNew(/*hyperplanes*/closeHyps, pyramids, pyr);


            pyramids.Add(pyr); //Add new pyramid to list
            //Add new pyramid's hyperplanes to list
            for (int j = 0; j < pyr.hyperplanes.Length; j++)
            {
                hyperplanes.Add(pyr.hyperplanes[j]);
            }

            //intersections.AddRange(newIntersections);
             //intersections = ih.PruneIntersections(intersections, pyramids); //Remove intersection points that lie below any of the pyramids
            newIntersections = ih.PruneIntersections(newIntersections, pyramids);
            intersections = ih.PruneIntersections(intersections, new List<Pyramid> { pyr });
            intersections.AddRange(newIntersections);


            //Remove past intersection dots
            Destroy(pointContainer);
            pointContainer = new GameObject();
            pointContainer.name = "Points";

            //Draw the intersections
            for (int i = 0; i < intersections.Count; i++)
            {
                dg.GraphPoint(pointContainer, intersections[i], 0.03f, Color.red);
            }

            dg.GraphPoint(pointContainer, Vector3.zero, 0.015f, Color.blue); //Solution

            print("Lower Bound: " + lowerBound);
            print("Upper Bound: " + upperBound);

            pyramidCount = pyramids.Count;
            hyperplaneCount = hyperplanes.Count;
        }

    }

    int MinSect(List<Vector3> lst)
    {
        float minVal = lst[0].y;
        int minPos = 0;

        for(int i = 1; i < lst.Count; i++)
        {
            if(lst[i].y < minVal)
            {
                minVal = lst[i].y;
                minPos = i;
            }
        }

        return minPos;
    }

    List<Hyperplane> GetAdjHyperplanes(List<Pyramid> pyrs, Pyramid pyr, int adjSize)
    {
        List<Hyperplane> adj = new List<Hyperplane>();

        for(int i = 0; i < pyrs.Count; i++)
        {
            pyrs[i].dist = Mathf.Max(Mathf.Abs(pyrs[i].peak.x - pyr.peak.x), Mathf.Abs(pyrs[i].peak.z - pyr.peak.z));
                //Vector2.Distance(new Vector2(pyrs[i].peak.x, pyrs[i].peak.z), new Vector2(pyr.peak.x, pyr.peak.z));
        }

        pyrs = pyrs.OrderBy(x => x.dist).ToList();
        for(int i = 0; i < Mathf.Min(adjSize, pyrs.Count); i++)
        {
            for(int j = 0; j < 4; j++)
            {
                adj.Add(pyrs[i].hyperplanes[j]);
            }           
        }

        return adj;
    }

    public float Objective(float x, float z)
    {
        return Objectives.QuadSine(x, z);
    }
}
