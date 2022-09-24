using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using mattatz.Triangulation2DSystem;
using DelaunatorSharp;
using DelaunatorSharp.Unity.Extensions;

public class PolygonTest : MonoBehaviour
{
    //mesh properties
    Vector3[] polygonPoints;
    int[] polygonTriangles;

    //polygon properties
    public int polygonSides;
    public float polygonRadius;
    public float centerRadius;

    public int modes = 4;
    public float[] modeCoeffs = new float[4];
    public float[] modeSinCoeffs = new float[4];
    public float area;
    public float realArea;
    public float offset = 0f;

    public LineRenderer line;
    public int lineRes = 10;

    public GameObject dotPrefab;
    List<GameObject> dots = new List<GameObject>();
    // float dotDist = 0.5f;
    public int dotCount = 10;
    public Transform dotParent;

    public Slider[] CosSliders;
    public Slider[] SinSliders;

    public Text[] cosTexts;
    public Text[] sinTexts;

    Mesh mesh;
    void Start()
    {
        // modeCoeffs = new float[modes];
        // modeSinCoeffs = new float[modes];
        // for (int i = 0; i < modes; i++)
        // {
        //     // modeCoeffs[i] = 0f;
        //     // modeSinCoeffs[i] = 0f;
        // }
        modeCoeffs[0] = polygonRadius;
        
        // DrawFilled(polygonSides,polygonRadius);
    }
    
    void DrawLine()
    {
        List<Vector3> points = new List<Vector3>();
        float circumferenceProgressPerStep = (float)1/polygonSides;
        float TAU = 2*Mathf.PI;
        float radianProgressPerStep = circumferenceProgressPerStep*TAU;
        float scale;
        float sqrdSum = 0f;
        for (int i = 0; i < polygonSides; i++)
        {
            float currentRadian = radianProgressPerStep*i;
            float newRadius = 0f;
            for (int j = 0; j < modes; j++)
            {
                newRadius += modeCoeffs[j]*Mathf.Cos(j*(currentRadian+offset)) + modeSinCoeffs[j]*Mathf.Sin(j*(currentRadian+offset));
            }
            sqrdSum += newRadius * newRadius;
        }
        scale = Mathf.Sqrt((2*area)/(sqrdSum * polygonSides * Mathf.Sin(TAU/polygonSides)));

        realArea = 0f;
        for(int i = 0; i < polygonSides; i++)
        {
            float currentRadian = radianProgressPerStep*i;
            Vector3 newPoint = new Vector3(Mathf.Cos(currentRadian), Mathf.Sin(currentRadian),0);
            float newRadius = 0f;
            for (int j = 0; j < modes; j++)
            {
                newRadius += modeCoeffs[j]*Mathf.Cos(j*(currentRadian+offset)) + modeSinCoeffs[j]*Mathf.Sin(j*(currentRadian+offset));
            }
            newRadius *= scale;
            newPoint *= newRadius;
            realArea += newRadius*newRadius;
            points.Add(newPoint);
        }
        realArea *= (polygonSides * Mathf.Sin(TAU/polygonSides))/2f;

        float perimeterLength = 2f * Mathf.PI * modeCoeffs[0] * scale;
        float segmentLength = perimeterLength/dotCount;
        float theta = 0f;
        float dtheta,r;
        int spawnedDots = 0;
        Vector3 spawnPos;
        while (theta <= 2*Mathf.PI)
        {
            r = 0f;
            for (int i = 0; i < modes; i++)
            {
                r += modeCoeffs[i]*Mathf.Cos(i*theta) + modeSinCoeffs[i]*Mathf.Sin(i*theta);
            }
            r *= scale;
            if(r <= Mathf.Epsilon)
            {
                theta += 0.1f * Mathf.PI;
                continue;
            }
            spawnPos = new Vector3(Mathf.Cos(theta), Mathf.Sin(theta),0);
            spawnPos *= r;
            if(dots.Count < spawnedDots + 1)
            {
                dots.Add(
                    Instantiate(dotPrefab,spawnPos,Quaternion.identity,dotParent)
                );
            }
            else
            {
                dots[spawnedDots].transform.position = spawnPos;
            }
            spawnedDots++;

            dtheta = segmentLength/r;
            theta += dtheta;
        }
        if(spawnedDots<dots.Count)
        {
            for (int i = 0; i < dots.Count - spawnedDots; i++)
            {
                Destroy(dots[dots.Count-1]);
                dots.RemoveAt(dots.Count-1);
            }
        }
        points.Add(points[0]);
        line.positionCount = points.Count;
        line.SetPositions(points.ToArray());
    }

    
    void Update()
    {
        // DrawLine();
        // Vector3[] linePos = new Vector3[lineRes+1];
        // for (int i = 0; i <= lineRes; i++)
        // {
        //     float x = (i*16f)/(float)lineRes;
        //     float y = 0f;
        //     float currentRadian = (2f*Mathf.PI * x) / 16f;
        //     for (int j = 0; j < modes; j++)
        //     {
        //         y += modeCoeffs[j]*Mathf.Cos(j*(currentRadian)) + modeSinCoeffs[j]*Mathf.Sin(j*(currentRadian));
        //     }
        //     linePos[i] = new Vector3(x,y,0);
        // }
        // line.positionCount = lineRes;
        // line.SetPositions(linePos);
        
        // DrawFilled(polygonSides,polygonRadius);
    }

    public void OnSliderChange()
    {
        for (int i = 0; i < 3; i++)
        {
            modeCoeffs[i+1] = CosSliders[i].value;
            cosTexts[i].text = "a" + (i+1).ToString() + " " + "=" + " " + modeCoeffs[i+1].ToString("0.000");
        }
        for (int i = 0; i < 3; i++)
        {
            modeSinCoeffs[i+1] = SinSliders[i].value;
            sinTexts[i].text = "b" + (i+1).ToString() + " " + "=" + " " + modeSinCoeffs[i+1].ToString("0.000");
        }
        DrawLine();
    }

    void UpdateVertices()
    {
        // DrawLine();
        // mesh.vertices = GetCircumferencePointsVector3(polygonSides,polygonRadius).ToArray();
        // GetComponent<MeshFilter>().sharedMesh = mesh;
    }

    // void DrawFilled(int sides, float radius)
    // {
    //     List<IPoint> polygonPointsList = GetCircumferencePointsIPoint(sides,radius);
    //     delaunator = new Delaunator(polygonPointsList.ToArray());
    //     mesh = new Mesh
    //     {
    //         vertices = delaunator.Points.ToVectors3(),
    //         triangles = delaunator.Triangles
    //     };

    //     mesh.RecalculateNormals();
    //     mesh.RecalculateBounds();
    //     GetComponent<MeshFilter>().sharedMesh = mesh;
    //     // List<Vector2> polygonPointsList = GetCircumferencePoints(sides,radius);
    //     // // polygonPoints = polygonPointsList.ToArray();
    //     // // polygonTriangles = DrawFilledTriangles(polygonPoints);
    //     // Polygon2D polygon = Polygon2D.Contour(polygonPointsList.ToArray());
    //     // // delaunator = new Delaunator(points.ToArray());

    //     // // construct Triangulation2D with Polygon2D and threshold angle (18f ~ 27f recommended)
    //     // Triangulation2D triangulation = new Triangulation2D(polygon, 22.5f);

    //     // // build a mesh from triangles in a Triangulation2D instance
    //     // mesh = triangulation.Build();
    //     // GetComponent<MeshFilter>().sharedMesh = mesh;
    //     // // mesh.vertices = polygonPoints;
    //     // // mesh.triangles = polygonTriangles;
    // }

    List<Vector2> GetCircumferencePoints(int sides, float radius)   
    {
        List<Vector2> points = new List<Vector2>();
        float circumferenceProgressPerStep = (float)1/sides;
        float TAU = 2*Mathf.PI;
        float radianProgressPerStep = circumferenceProgressPerStep*TAU;
        float scale;
        float sqrdSum = 0f;
        for (int i = 0; i < sides; i++)
        {
            float currentRadian = radianProgressPerStep*i;
            float newRadius = 0f;
            for (int j = 0; j < modes; j++)
            {
                newRadius += modeCoeffs[j]*Mathf.Cos(j*(currentRadian+offset)) + modeSinCoeffs[j]*Mathf.Sin(j*(currentRadian+offset));
            }
            sqrdSum += newRadius * newRadius;
        }
        scale = Mathf.Sqrt((2*area)/(sqrdSum * sides * Mathf.Sin(TAU/sides)));

        realArea = 0f;
        for(int i = 0; i < sides; i++)
        {
            float currentRadian = radianProgressPerStep*i;
            Vector2 newPoint = new Vector2(Mathf.Cos(currentRadian), Mathf.Sin(currentRadian));
            float newRadius = 0f;
            for (int j = 0; j < modes; j++)
            {
                newRadius += modeCoeffs[j]*Mathf.Cos(j*(currentRadian+offset)) + modeSinCoeffs[j]*Mathf.Sin(j*(currentRadian+offset));
            }
            newRadius *= scale;
            newPoint *= newRadius;
            realArea += newRadius*newRadius;
            points.Add(newPoint);
        }
        realArea *= (sides * Mathf.Sin(TAU/sides))/2f;

        float perimeterLength = 2f * Mathf.PI * modeCoeffs[0] * scale;
        float segmentLength = perimeterLength/dotCount;
        float theta = 0f;
        float dtheta,r;
        int spawnedDots = 0;
        Vector3 spawnPos;
        while (theta <= 2*Mathf.PI)
        {
            r = 0f;
            for (int i = 0; i < modes; i++)
            {
                r += modeCoeffs[i]*Mathf.Cos(i*theta) + modeSinCoeffs[i]*Mathf.Sin(i*theta);
            }
            r *= scale;
            spawnPos = new Vector3(Mathf.Cos(theta), Mathf.Sin(theta),0);
            spawnPos *= r;
            if(dots.Count < spawnedDots + 1)
            {
                dots.Add(
                    Instantiate(dotPrefab,spawnPos,Quaternion.identity,dotParent)
                );
            }
            else
            {
                dots[spawnedDots].transform.position = spawnPos;
            }
            spawnedDots++;

            dtheta = segmentLength/r;
            theta += dtheta;
        }
        if(spawnedDots<dots.Count)
        {
            for (int i = 0; i < dots.Count - spawnedDots; i++)
            {
                Destroy(dots[dots.Count-1]);
                dots.RemoveAt(dots.Count-1);
            }
        }
        return points;
    }

    List<Vector3> GetCircumferencePointsVector3(int sides, float radius)   
    {
        List<Vector3> points = new List<Vector3>();
        float circumferenceProgressPerStep = (float)1/sides;
        float TAU = 2*Mathf.PI;
        float radianProgressPerStep = circumferenceProgressPerStep*TAU;
        float scale;
        float sqrdSum = 0f;
        for (int i = 0; i < sides; i++)
        {
            float currentRadian = radianProgressPerStep*i;
            float newRadius = 0f;
            for (int j = 0; j < modes; j++)
            {
                newRadius += modeCoeffs[j]*Mathf.Cos(j*(currentRadian+offset)) + modeSinCoeffs[j]*Mathf.Sin(j*(currentRadian+offset));
            }
            sqrdSum += newRadius * newRadius;
        }
        scale = Mathf.Sqrt((2*area)/(sqrdSum * sides * Mathf.Sin(TAU/sides)));

        realArea = 0f;
        for(int i = 0; i < sides; i++)
        {
            float currentRadian = radianProgressPerStep*i;
            Vector3 newPoint = new Vector3(Mathf.Cos(currentRadian), Mathf.Sin(currentRadian),0);
            float newRadius = 0f;
            for (int j = 0; j < modes; j++)
            {
                newRadius += modeCoeffs[j]*Mathf.Cos(j*(currentRadian+offset)) + modeSinCoeffs[j]*Mathf.Sin(j*(currentRadian+offset));
            }
            newRadius *= scale;
            newPoint *= newRadius;
            realArea += newRadius*newRadius;
            points.Add(newPoint);
        }
        realArea *= (sides * Mathf.Sin(TAU/sides))/2f;

        float perimeterLength = 2f * Mathf.PI * modeCoeffs[0] * scale;
        float segmentLength = perimeterLength/dotCount;
        float theta = 0f;
        float dtheta,r;
        int spawnedDots = 0;
        Vector3 spawnPos;
        while (theta <= 2*Mathf.PI)
        {
            r = 0f;
            for (int i = 0; i < modes; i++)
            {
                r += modeCoeffs[i]*Mathf.Cos(i*theta) + modeSinCoeffs[i]*Mathf.Sin(i*theta);
            }
            r *= scale;
            if(r <= Mathf.Epsilon)
            {
                theta += 0.1f * Mathf.PI;
                continue;
            }
            spawnPos = new Vector3(Mathf.Cos(theta), Mathf.Sin(theta),0);
            spawnPos *= r;
            if(dots.Count < spawnedDots + 1)
            {
                dots.Add(
                    Instantiate(dotPrefab,spawnPos,Quaternion.identity,dotParent)
                );
            }
            else
            {
                dots[spawnedDots].transform.position = spawnPos;
            }
            spawnedDots++;

            dtheta = segmentLength/r;
            theta += dtheta;
        }
        if(spawnedDots<dots.Count)
        {
            for (int i = 0; i < dots.Count - spawnedDots; i++)
            {
                Destroy(dots[dots.Count-1]);
                dots.RemoveAt(dots.Count-1);
            }
        }
        return points;
    }

    List<IPoint> GetCircumferencePointsIPoint(int sides, float radius)   
    {
        List<IPoint> points = new List<IPoint>();
        float circumferenceProgressPerStep = (float)1/sides;
        float TAU = 2*Mathf.PI;
        float radianProgressPerStep = circumferenceProgressPerStep*TAU;
        float scale;
        float sqrdSum = 0f;
        for (int i = 0; i < sides; i++)
        {
            float currentRadian = radianProgressPerStep*i;
            float newRadius = 0f;
            for (int j = 0; j < modes; j++)
            {
                newRadius += modeCoeffs[j]*Mathf.Cos(j*(currentRadian+offset)) + modeSinCoeffs[j]*Mathf.Sin(j*(currentRadian+offset));
            }
            sqrdSum += newRadius * newRadius;
        }
        scale = Mathf.Sqrt((2*area)/(sqrdSum * sides * Mathf.Sin(TAU/sides)));

        realArea = 0f;
        for(int i = 0; i < sides; i++)
        {
            float currentRadian = radianProgressPerStep*i;
            float newRadius = 0f;
            for (int j = 0; j < modes; j++)
            {
                newRadius += modeCoeffs[j]*Mathf.Cos(j*(currentRadian+offset)) + modeSinCoeffs[j]*Mathf.Sin(j*(currentRadian+offset));
            }
            newRadius *= scale;
            realArea += newRadius*newRadius;
            points.Add(new Point(Mathf.Cos(currentRadian) * newRadius, Mathf.Sin(currentRadian) * newRadius));
        }
        realArea *= (sides * Mathf.Sin(TAU/sides))/2f;

        float perimeterLength = 2f * Mathf.PI * modeCoeffs[0] * scale;
        float segmentLength = perimeterLength/dotCount;
        float theta = 0f;
        float dtheta,r;
        int spawnedDots = 0;
        Vector3 spawnPos;
        while (theta <= 2*Mathf.PI)
        {
            r = 0f;
            for (int i = 0; i < modes; i++)
            {
                r += modeCoeffs[i]*Mathf.Cos(i*theta) + modeSinCoeffs[i]*Mathf.Sin(i*theta);
            }
            r *= scale;
            spawnPos = new Vector3(Mathf.Cos(theta), Mathf.Sin(theta),0);
            spawnPos *= r;
            if(dots.Count < spawnedDots + 1)
            {
                dots.Add(
                    Instantiate(dotPrefab,spawnPos,Quaternion.identity,dotParent)
                );
            }
            else
            {
                dots[spawnedDots].transform.position = spawnPos;
            }
            spawnedDots++;

            dtheta = segmentLength/r;
            theta += dtheta;
        }
        if(spawnedDots<dots.Count)
        {
            for (int i = 0; i < dots.Count - spawnedDots; i++)
            {
                Destroy(dots[dots.Count-1]);
                dots.RemoveAt(dots.Count-1);
            }
        }
        return points;
    }
    
    int[] DrawFilledTriangles(Vector3[] points)
    {   
        int triangleAmount = points.Length - 2;
        List<int> newTriangles = new List<int>();
        for(int i = 0; i<triangleAmount; i++)
        {
            newTriangles.Add(0);
            newTriangles.Add(i+1);
            newTriangles.Add(i+2);
        }
        return newTriangles.ToArray();
    }
}



