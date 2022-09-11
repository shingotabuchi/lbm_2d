using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class D2Q4Main : MonoBehaviour
{
    public Image plotImage;
    public int resInt;
    [Range(0.25f, 2.0f)]
    public float alpha;
    public BoundaryType[] wallboundaries = new BoundaryType[]{BoundaryType.Constant,BoundaryType.Constant,BoundaryType.Constant,BoundaryType.Constant};
    [Range(0.0f, 1.0f)]
    public float wallTemp1 = 1f;
    [Range(0.0f, 1.0f)]
    public float wallTemp2 = 0f;
    [Range(0.0f, 1.0f)]
    public float wallTemp3 = 0f;
    [Range(0.0f, 1.0f)]
    public float wallTemp4 = 0f;

    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    float[,] f1,f2,f3,f4,rho;
    float dx = 1f;
    float dt = 1f;
    float csq,omega,feq;
    public int loopCount = 1;
    void Start()
    {
        plotTexture = new Texture2D(resInt,resInt);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,resInt,resInt),Vector2.zero);
        csq=dx*dx/(dt*dt);
        f1 = new float[resInt,resInt];
        f2 = new float[resInt,resInt];
        f3 = new float[resInt,resInt];
        f4 = new float[resInt,resInt];
        rho = new float[resInt,resInt];
        for (int i = 0; i < resInt; i++)
        {
            for (int j = 0; j < resInt; j++)
            {
                rho[i,j] = 0f;
                f1[i,j] = 0.25f*rho[i,j];
                f2[i,j] = 0.25f*rho[i,j];
                f3[i,j] = 0.25f*rho[i,j];
                f4[i,j] = 0.25f*rho[i,j];
            }
        }
        UpdatePlot();
    }

    void LBMStep()
    {
        // Collision Process
        for (int i = 0; i < resInt; i++)
        {
            for (int j = 0; j < resInt; j++)
            {
                feq = 0.25f*rho[i,j];
                f1[i,j] = omega*feq + (1f - omega)*f1[i,j];
                f2[i,j] = omega*feq + (1f - omega)*f2[i,j];
                f3[i,j] = omega*feq + (1f - omega)*f3[i,j];
                f4[i,j] = omega*feq + (1f - omega)*f4[i,j];
            }
        }
        // Streaming Process
        for (int j = 0; j < resInt; j++)
        {
            for (int i = 1; i < resInt; i++)
            {
                f1[resInt - i,j] = f1[resInt - i - 1,j];
                f2[i-1,j] = f2[i,j];
            }
        }
        for (int i = 0; i < resInt; i++)
        {
            for (int j = 1; j < resInt; j++)
            {
                f3[i,resInt - j] = f3[i,resInt - j - 1];
                f4[i,j-1] = f4[i,j];
            }
        }
        // Boundary Conditions
        for (int j = 0; j < resInt; j++)
        {
            switch (wallboundaries[0])
            {
                case BoundaryType.Constant: 
                f1[0,j] = 0.5f*wallTemp1 - f2[0,j];
                f3[0,j] = 0.5f*wallTemp1 - f4[0,j];
                break;
                case BoundaryType.Adiabatic:
                f1[0,j] = f1[1,j];
                f2[0,j] = f2[1,j];
                f3[0,j] = f3[1,j];
                f4[0,j] = f4[1,j];
                break;
            }
            switch (wallboundaries[1])
            {
                case BoundaryType.Constant: 
                f2[resInt-1,j] = 0.5f*wallTemp2 - f1[resInt-1,j];
                f3[resInt-1,j] = 0.5f*wallTemp2 - f4[resInt-1,j];
                break;
                case BoundaryType.Adiabatic:
                f1[resInt-1,j] = f1[resInt-2,j];
                f2[resInt-1,j] = f2[resInt-2,j];
                f3[resInt-1,j] = f3[resInt-2,j];
                f4[resInt-1,j] = f4[resInt-2,j];
                break;
            }
        }
        for (int i = 1; i < resInt; i++)
        {
            switch (wallboundaries[2])
            {
                case BoundaryType.Constant: 
                f1[i,resInt-1] = 0.5f*wallTemp3 - f2[i,resInt-1];
                f4[i,resInt-1] = 0.5f*wallTemp3 - f3[i,resInt-1];
                break;
                case BoundaryType.Adiabatic:
                f1[i,resInt-1] = f1[i,resInt-2];
                f2[i,resInt-1] = f2[i,resInt-2];
                f3[i,resInt-1] = f3[i,resInt-2];
                f4[i,resInt-1] = f4[i,resInt-2];
                break;
            }
            switch (wallboundaries[3])
            {
                case BoundaryType.Constant: 
                f1[i,0] = 0.5f*wallTemp4 - f2[i,0];
                f3[i,0] = 0.5f*wallTemp4 - f4[i,0];
                break;
                case BoundaryType.Adiabatic:
                f1[i,0] = f1[i,1];
                f2[i,0] = f2[i,1];
                f3[i,0] = f3[i,1];
                f4[i,0] = f4[i,1];
                break;
            }
        }
        // Get Rho
        for (int i = 0; i < resInt; i++)
        {
            for (int j = 0; j < resInt; j++)
            {
                rho[i,j]=f1[i,j]+f2[i,j]+f3[i,j]+f4[i,j];
            }
        }
    }
    
    void Update()
    {
        omega = 1f/(2f*alpha/(dt*csq)+0.5f);
        for (int i = 0; i < loopCount; i++)
        {
            LBMStep();
        }
        UpdatePlot();
    }

    void UpdatePlot()
    {
        for (int i = 0; i < plotPixels.Length; i++)
        {
            plotPixels[i] = colorHeatMap.GetColorForValue(rho[i%resInt,i/resInt],1f);
        }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
    }
}
