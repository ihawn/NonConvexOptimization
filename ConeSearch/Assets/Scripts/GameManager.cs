using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

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
    public int maxIterations = 5;
    public int pyramidCount = 0;

    public float lowerBound = -100000, upperBound = 100000;

    void Start()
    {
        UnityEditor.SceneView.FocusWindowIfItsOpen(typeof(UnityEditor.SceneView));

        pointContainer = new GameObject();
        hyperplanes = new List<Hyperplane>();


        Pyramid pyr1 = ph.GeneratePyramid(new Vector3(-0.8f, 0, -0.6f), L, 0);
        Pyramid pyr2 = ph.GeneratePyramid(new Vector3(-0.8f, 0, 0.95f), L, 1);
        Pyramid pyr3 = ph.GeneratePyramid(new Vector3(0.9f, 0, -0.6f), L, 2);
        Pyramid pyr4 = ph.GeneratePyramid(new Vector3(0.9f, 0, 0.95f), L, 3);

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
        for(int i = 0; i < pyramids.Count; i++)
        {
            for(int j = 0; j < pyramids[i].hyperplanes.Length; j++)
            {
                hyperplanes.Add(pyramids[i].hyperplanes[j]);
            }
        }
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


            //Combinations of pyramids intersecting
            /*combos = ph.CombinePyramids(pyramids);
            intersections = new List<Vector3>();
            for (int i = 0; i < combos.Count; i++)
            {
                intersections.AddRange(ih.IntersectPyramids(combos[i][0], combos[i][1], combos[i][2])); //Combine the lists
            }*/

            intersections.AddRange(ih.IntersectNew(hyperplanes, pyramids, intersections, pyr));
            intersections = intersections.Distinct().ToList();
            pyramids.Add(pyr); //Add new pyramid to list

            //Add new pyramid's hyperplanes to list
            for (int j = 0; j < pyr.hyperplanes.Length; j++)
            {
                hyperplanes.Add(pyr.hyperplanes[j]);
            }

            ih.PruneIntersections(intersections, pyramids); //Remove intersection points that lie below any of the pyramids

            //Remove past intersections
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

    public float Objective(float x, float z)
    {
        return x * x + z * z;// Objectives.Rastrigin(x, z);
    }
}
