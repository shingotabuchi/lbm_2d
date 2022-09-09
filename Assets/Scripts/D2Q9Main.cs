using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class D2Q9Main : MonoBehaviour
{
    public Image plotImage;
    public int resInt;
    [Range(0.25f, 2.0f)]
    public float alpha = 1f;
    public BoundaryType[] wallboundaries = new BoundaryType[]{BoundaryType.Constant,BoundaryType.Constant,BoundaryType.Constant,BoundaryType.Constant};
    [Range(0.0f, 1.0f)]
    public float wallTemp1 = 1f;
    [Range(0.0f, 1.0f)]
    public float wallTemp2 = 0f;
    [Range(0.0f, 1.0f)]
    public float wallTemp3 = 0f;
    [Range(0.0f, 1.0f)]
    public float wallTemp4 = 0f;
    public float u = 0.1f;
    public float v = 0.4f;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    float[,] rho;
    float[] w = new float[9];
    float[] feq = new float[9];
    float[,,] f;

    float dx = 1f;
    float dt = 1f;
    float csq,omega,ck;
    void Start()
    {
        plotTexture = new Texture2D(resInt,resInt);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,resInt,resInt),Vector2.zero);
        
        csq=dx*dx/(dt*dt);
        ck = dx/dt;
        rho = new float[resInt,resInt];
        f = new float[9,resInt,resInt];
        for (int i = 0; i < 9; i++)
        {
            if(i==0) w[i] = 4f/9f;
            else if(i<5) w[i] = 1f/9f;
            else w[i] = 1f/36f;
        }
        float initialDensity = 0f;
        for (int i = 0; i < resInt; i++)
        {
            for (int j = 0; j < resInt; j++)
            {
                float sum = 0f;
                for (int k = 0; k < 9; k++)
                {
                    f[k,i,j] = w[k] * initialDensity;
                    if(i==0) f[k,i,j] = w[k]*wallTemp1;
                    sum += f[k,i,j];
                }
                rho[i,j] = sum;
            }
        }
        UpdatePlot();
    }
    
    void Update()
    {

        omega = 1f/(3f*alpha/(dt*csq)+0.5f);
        // Collision Process
        for (int i = 0; i < resInt; i++)
        {
            for (int j = 0; j < resInt; j++)
            {
                feq[0] = w[0] * rho[i,j];
                feq[1] = w[1] * rho[i,j] * (1f + 3.0f*u/ck);
                feq[2] = w[2] * rho[i,j] * (1f + 3.0f*v/ck);
                feq[3] = w[3] * rho[i,j] * (1f - 3.0f*u/ck);
                feq[4] = w[4] * rho[i,j] * (1f - 3.0f*v/ck);
                feq[5] = w[5] * rho[i,j] * (1f + 3.0f*(u+v)/ck);
                feq[6] = w[6] * rho[i,j] * (1f + 3.0f*(-u+v)/ck);
                feq[7] = w[7] * rho[i,j] * (1f + 3.0f*(-u-v)/ck);
                feq[8] = w[8] * rho[i,j] * (1f + 3.0f*(u-v)/ck);

                for (int k = 0; k < 9; k++)
                {
                    f[k,i,j] = (1 - omega) * f[k,i,j] + omega * feq[k];
                }
            }
        }
        // Streaming Process
        for (int j = resInt-1; j >= 0; j--)
        {
            for (int i = resInt-1; i >= 1; i--)
            {
                f[1,i,j] = f[1,i-1,j];
                if(j!=0) f[5,i,j] = f[5,i-1,j-1];
            }
        }
        for (int j = resInt-1; j >= 1; j--)
        {
            for (int i = 0; i < resInt; i++)
            {
                f[2,i,j] = f[2,i,j-1];
                if(i!=resInt-1)f[6,i,j] = f[6,i+1,j-1];
            }
        }
        for (int j = 0; j < resInt; j++)
        {
            for (int i = 0; i < resInt-1; i++)
            {
                f[3,i,j] = f[3,i+1,j];
                if(j!=resInt-1)f[7,i,j] = f[7,i+1,j+1];
            }
        }
        for (int j = 0; j < resInt-1; j++)
        {
            for (int i = resInt-1; i >= 0; i--)
            {
                f[4,i,j] = f[4,i,j+1];
                if(i!=0)f[8,i,j] = f[8,i-1,j+1];
            }
        }
        // Boundary Conditions
        for (int j = 0; j < resInt; j++)
        {
            switch (wallboundaries[0])
            {
                case BoundaryType.Constant: 
                f[1,0,j] = (w[1] + w[3])*wallTemp1 - f[3,0,j];
                f[5,0,j] = (w[5] + w[7])*wallTemp1 - f[7,0,j];
                f[8,0,j] = (w[8] + w[6])*wallTemp1 - f[6,0,j];
                f[2,0,j] = (w[2] + w[4])*wallTemp1 - f[4,0,j];//なかった
                f[0,0,j] = w[0]*wallTemp1;//なかった
                break;
                case BoundaryType.Adiabatic:
                for (int k = 0; k < 9; k++)
                {
                    f[k,0,j] = f[k,1,j];
                }
                break;
            }
            switch (wallboundaries[1])
            {
                case BoundaryType.Constant: 
                f[3,resInt-1,j] = (w[3] + w[1])*wallTemp2 - f[1,resInt-1,j];
                f[6,resInt-1,j] = (w[6] + w[8])*wallTemp2 - f[8,resInt-1,j];
                f[7,resInt-1,j] = (w[7] + w[5])*wallTemp2 - f[5,resInt-1,j];
                f[2,resInt-1,j] = (w[2] + w[4])*wallTemp2 - f[4,resInt-1,j];
                f[0,resInt-1,j] = w[0]*wallTemp2;
                break;
                case BoundaryType.Adiabatic:
                for (int k = 0; k < 9; k++)
                {
                    f[k,resInt-1,j] = f[k,resInt-2,j];
                }
                break;
            }
        }
        for (int i = 1; i < resInt; i++)
        {
            switch (wallboundaries[2])
            {
                case BoundaryType.Constant: 
                f[4,i,resInt-1] = (w[4] + w[2])*wallTemp3 - f[2,i,resInt-1];
                f[7,i,resInt-1] = (w[7] + w[5])*wallTemp3 - f[5,i,resInt-1];
                f[8,i,resInt-1] = (w[8] + w[6])*wallTemp3 - f[6,i,resInt-1];
                f[1,i,resInt-1] = (w[1] + w[3])*wallTemp3 - f[3,i,resInt-1];
                f[0,i,resInt-1] = w[0]*wallTemp3;
                break;
                case BoundaryType.Adiabatic:
                for (int k = 0; k < 9; k++)
                {
                    f[k,i,resInt-1] = f[k,i,resInt-2];
                }
                break;
            }
            switch (wallboundaries[3])
            {
                case BoundaryType.Constant: 
                f[2,i,0] = (w[2] + w[4])*wallTemp4 - f[4,i,0];
                f[5,i,0] = (w[5] + w[7])*wallTemp4 - f[7,i,0];
                f[6,i,0] = (w[6] + w[8])*wallTemp4 - f[8,i,0];
                f[1,i,0] = (w[1] + w[3])*wallTemp4 - f[3,i,0];
                f[0,i,0] = w[0]*wallTemp4;
                break;
                case BoundaryType.Adiabatic:
                for (int k = 0; k < 9; k++)
                {
                    f[k,i,0] = f[k,i,1];
                }
                break;
            }
        }
        // Get Rho
        for (int i = 0; i < resInt; i++)
        {
            for (int j = 0; j < resInt; j++)
            {
                float sum = 0f;
                for (int k = 0; k < 9; k++)
                {
                    sum += f[k,i,j];
                }
                if(sum>1f) sum = 1f; // 微妙に１より大きい時がある。原因は謎。
                rho[i,j]=sum;
            }
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
