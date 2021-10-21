using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;

public class Hyperplane
{
    public int parentID;
    public int direction; //0-3
    public Vector<double> coeff;

    public Hyperplane(int parentID, int direction, Vector<double> coeff)
    {
        this.parentID = parentID;
        this.direction = direction;
        this.coeff = coeff;
    }

    public Hyperplane(Hyperplane hyperplane)
    {
        this.parentID = hyperplane.parentID;
        this.direction = hyperplane.direction;
        this.coeff = hyperplane.coeff;
    }
}
