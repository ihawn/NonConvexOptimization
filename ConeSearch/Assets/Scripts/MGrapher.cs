using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MGrapher : MonoBehaviour
{
    public PyramidHandler pHandler;
    public GameManager gm;

    public void GeneratePyramidMesh(Pyramid p)
    {
        Vector3 peak = p.peak;
        peak.y = gm.Objective(peak.x, peak.z);
        float L = p.L;

        Vector3[] pts = new Vector3[]
        {
            new Vector3(peak.x + 1, peak.y - L, peak.z + 1),
            new Vector3(peak.x - 1, peak.y - L, peak.z + 1),
            new Vector3(peak.x - 1, peak.y - L, peak.z - 1),
            new Vector3(peak.x + 1, peak.y - L, peak.z - 1),
            peak
        };


        //Displaying the mesh
        GameObject g = new GameObject();
        g.name = "Pyramid";
        MeshRenderer mr = g.AddComponent<MeshRenderer>();
        MeshFilter mf = g.AddComponent<MeshFilter>();
        mr.sharedMaterial = new Material(Shader.Find("Standard"));
        Mesh mesh = new Mesh();

        Vector3[] verts = new Vector3[12]
        {
            pts[0], pts[4], pts[1],
            pts[1], pts[4], pts[2],
            pts[2], pts[4], pts[3],
            pts[3], pts[4], pts[0]
        };

        int[] tris = new int[12] {
            0, 1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 11
        };

        mesh.vertices = verts;
        mesh.triangles = tris;
        mf.mesh = mesh;
        mesh.RecalculateNormals();
        mesh.Optimize();

        //Random color for each pyramid
        Color col = new Color(Random.Range(0f, 1f), Random.Range(0f, 1f), Random.Range(0f, 1f));
        mr.material.SetColor("_Color", col);
    }
}
