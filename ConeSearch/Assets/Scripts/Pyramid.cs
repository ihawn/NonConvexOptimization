using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Pyramid
{
    public int id;
    public Vector3 peak;
    public float L;
    public float dist; //distance to another pyramid in question
    public Hyperplane[] hyperplanes;

    public Pyramid(int id, Vector3 peak, float L, float dist, Hyperplane[] hyperplanes)
    {
        this.id = id;
        this.peak = peak;
        this.L = L;
        this.dist = dist;
        this.hyperplanes = hyperplanes;
    }

    public Pyramid(Pyramid pyramid)
    {
        this.id = pyramid.id;
        this.peak = pyramid.peak;
        this.L = pyramid.L;
        this.dist = pyramid.dist;
        this.hyperplanes = pyramid.hyperplanes;
    }
}
