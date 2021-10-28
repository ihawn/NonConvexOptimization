using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class Objectives
{
    public static float StyblinskiTang(float x, float z)
    {
        return (Mathf.Pow(x, 4) - 16 * Mathf.Pow(x, 2) + 5 * x + Mathf.Pow(z, 4) - 16 * Mathf.Pow(z, 2) + 5 * z) / 2f;
    }

    public static float Rastrigin(float x, float z)
    {
        return 20 + Mathf.Pow(x, 2) - 10 * Mathf.Cos(2 * Mathf.PI * x) + Mathf.Pow(z, 2) - 10 * Mathf.Cos(2 * Mathf.PI * z);
    }

    public static float QuadSine(float x, float y)
    {
        return 0.05f * (x * x + y * y) + 0.5f * Mathf.Sin((x - 1) * y);
    }
}
