using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class ChannelFlow : MonoBehaviour
{
    public Image plotImage;
    public int resInt;
    int resIntX,resIntY;
    [Range(0.005f, 0.1f)]
    public float alpha = 0.01f;
    public float aspectRatio = 25f;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    float[,] rho;
    float[] w = new float[9];
    float[] cx = new float[9]{0, 1, 0,-1, 0, 1,-1,-1, 1};
    float[] cy = new float[9]{0, 0, 1, 0,-1, 1, 1,-1,-1};
    float[,,] feq;
    float[,,] f,ftmp;
    float[,] u;
    float[,] v;
    float[,] speed;
    float maxSpeed;
    [Range(0.01f, 0.2f)]
    public float u0 = 0.1f;
    public VectorField vectorField;
    public float rho0 = 5f;
    public HeatMapMode mode = HeatMapMode.Speed;
    float maxrho;
    float omega,Re;
    void Start()
    {
        resIntX = (int)(resInt*aspectRatio);
        resIntY = resInt;
        plotTexture = new Texture2D(resIntX,resIntY);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,resIntX,resIntY),Vector2.zero);
        
        Re = u0*resInt/alpha;
        print(Re);
        rho = new float[resIntX,resIntY];
        u = new float[resIntX,resIntY];
        v = new float[resIntX,resIntY];
        speed = new float[resIntX,resIntY];
        feq = new float[9,resIntX,resIntY];
        f = new float[9,resIntX,resIntY];
        ftmp = new float[9,resIntX,resIntY];
        for (int i = 0; i < 9; i++)
        {
            if(i==0) w[i] = 4f/9f;
            else if(i<5) w[i] = 1f/9f;
            else w[i] = 1f/36f;
        }
        for (int i = 0; i < resIntX; i++)
        {
            for (int j = 0; j < resIntY; j++)
            {
                rho[i,j] = rho0;
                u[i,j] = 0f;
                v[i,j] = 0f;
                speed[i,j] = 0f;
                if(i == 0)
                {
                    u[i,j] = u0;
                    v[i,j] = 0f;
                    speed[i,j] = u0;
                }
            }
        }
        maxrho = rho0;
        maxSpeed = u0;
        UpdatePlot();
    }
    
    void Update()
    {
        omega = 1f/(3f*alpha+0.5f);
        Collisions();
        Streaming();
        SfBound();
        RhoUV();
        vectorField.UpdateVectors(u,v,resInt,maxSpeed);
        UpdatePlot();
    }

    void Collisions()
    {
        for (int i = 0; i < resIntX; i++)
        {
            for (int j = 0; j < resInt; j++)
            {
                float t1 = u[i,j]*u[i,j] + v[i,j]*v[i,j];
                for (int k = 0; k < 9; k++)
                {
                    float t2 = u[i,j]*cx[k] + v[i,j]*cy[k];
                    feq[k,i,j] = rho[i,j]*w[k]*(1.0f+3.0f*t2+4.50f*t2*t2-1.50f*t1);
                    f[k,i,j]=omega*feq[k,i,j]+(1f-omega)*f[k,i,j];
                }
            }
        }
    }

    void Streaming()
    {
        ftmp = (float[,,])(f.Clone());
        for(int i = 0; i < resIntX; i++)
        { 
            for(int j = 0; j < resInt; j++)
            { 
                for(int k = 0; k <= 8; k++)
                {
                    //periodic boundary condition
                    int im = i + (int)cx[k]; 
                    int jm = j + (int)cy[k];
                    if(
                        (jm != resInt && jm !=-1)
                        &&
                        (im != resIntX && im !=-1)
                    ) f[k,im,jm] = ftmp[k,i,j];
                } 
            }
        }
        // for (int j = resInt-1; j >= 0; j--)
        // {
        //     for (int i = resIntX-1; i >= 1; i--)
        //     {
        //         f[1,i,j] = f[1,i-1,j];
        //         if(j!=0) f[5,i,j] = f[5,i-1,j-1];
        //     }
        // }
        // for (int j = resInt-1; j >= 1; j--)
        // {
        //     for (int i = 0; i < resIntX; i++)
        //     {
        //         f[2,i,j] = f[2,i,j-1];
        //         if(i!=resIntX-1)f[6,i,j] = f[6,i+1,j-1];
        //     }
        // }
        // for (int j = 0; j < resInt; j++)
        // {
        //     for (int i = 0; i < resIntX-1; i++)
        //     {
        //         f[3,i,j] = f[3,i+1,j];
        //         if(j!=resInt-1)f[7,i,j] = f[7,i+1,j+1];
        //     }
        // }
        // for (int j = 0; j < resInt-1; j++)
        // {
        //     for (int i = resIntX-1; i >= 0; i--)
        //     {
        //         f[4,i,j] = f[4,i,j+1];
        //         if(i!=0)f[8,i,j] = f[8,i-1,j+1];
        //     }
        // }
    }

    void SfBound()
    {
        for (int j = 0; j < resInt; j++)
        {
            // Left
            float localRho = (f[0,0,j] + f[2,0,j] + f[4,0,j] + 2f*(f[3,0,j] + f[6,0,j] + f[7,0,j]))/(1f - u0);
            f[1,0,j] = f[3,0,j] + 2f*localRho * u0/3f;
            f[5,0,j] = f[7,0,j] + localRho * u0/6;
            f[8,0,j] = f[6,0,j] + localRho * u0/6;
            // Right
            f[1,resIntX-1,j] = 2f*f[1,resIntX-2,j] - f[1,resIntX-3,j];
            f[5,resIntX-1,j] = 2f*f[5,resIntX-2,j] - f[5,resIntX-3,j];
            f[8,resIntX-1,j] = 2f*f[8,resIntX-2,j] - f[8,resIntX-3,j];
        }
        for (int i = 0; i < resIntX; i++)
        {
            // Bottom
            f[2,i,0] = f[4,i,0];
            f[5,i,0] = f[7,i,0];
            f[6,i,0] = f[8,i,0];
            // Top
            f[4,i,resInt-1]=f[2,i,resInt-1];
            f[8,i,resInt-1]=f[6,i,resInt-1];
            f[7,i,resInt-1]=f[5,i,resInt-1];
        }
    }

    void RhoUV()
    {
        maxSpeed = 0f;
        maxrho = 0f;
        for (int i = 0; i < resIntX; i++)
        {
            for (int j = 0; j < resInt; j++)
            {

                if(i == 0)
                {
                    rho[i,j] = f[0,i,j] + f[1,i,j] + f[3,i,j] + 2f*(f[2,i,j] + f[6,i,j] + f[5,i,j]);
                }
                else
                {
                    float sum = 0f;
                    for (int k = 0; k < 9; k++)
                    {
                        sum += f[k,i,j];
                    }
                    rho[i,j] = sum;
                }
                maxrho = Mathf.Max(maxrho,rho[i,j]);

                if(i == 0)
                {
                    u[i,j] = u0;
                    v[i,j] = 0f;
                    speed[i,j] = u0;
                }
                else
                {
                    float usum = 0f;
                    float vsum = 0f;
                    for (int k = 0; k < 9; k++)
                    {
                        usum += f[k,i,j] * cx[k];
                        vsum += f[k,i,j] * cy[k];
                    }
                    u[i,j] = usum/rho[i,j];
                    v[i,j] = vsum/rho[i,j];
                    if(i== resIntX -1) v[i,j] = 0f;
                    speed[i,j] = Mathf.Sqrt(u[i,j]*u[i,j] + v[i,j]*v[i,j]);
                }
                maxSpeed = Mathf.Max(maxSpeed,speed[i,j]);
            }
        }
    }

    void UpdatePlot()
    {
        for (int i = 0; i < plotPixels.Length; i++)
        {
            if(mode == HeatMapMode.Speed)
            {
                plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%resIntX,i/resIntX],maxSpeed);
            }
            else if(mode == HeatMapMode.Density)
            {
                plotPixels[i] = colorHeatMap.GetColorForValue(rho[i%resIntX,i/resIntX],maxrho);
            }
        }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
    }
}
