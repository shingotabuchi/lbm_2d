using System.Collections;
using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LBMCylinder : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public VectorField vectorField;
    public bool vectorFieldOn;
    public HeatMapMode mode = HeatMapMode.Speed;
    public bool fixHeatMapMax;
    public bool fixHeatMapMin;
    public float fixedMaxSpeed = 0.6f;
    public float fixedMaxRho = 1.01f;
    public float fixedMinSpeed = 0f;
    public float fixedMinRho = 0f;
    public int DIM_X = 51;
    public int DIM_Y = 51;

    public float maxRho,minRho;
    public float maxSpeed,minSpeed;

    float[] cx = new float[9]{0, 1,    0,   -1,    0,     1,    -1,    -1,     1};
    float[] cy = new float[9]{0, 0,    1,    0,   -1,     1,     1,    -1,    -1};
    float[] w = new float[9]{4f/9f,1f/9f,1f/9f,1f/9f,1f/9f,1f/36f,1f/36f,1f/36f,1f/36f};
    float[,] rho, u, v, speed, fx ,fy;
    float[,,] f, f0, ftmp;
    public float u0 = 0.01f;
    public float tau = 0.6f;
    float nu;
    public int loopCount = 1;

    RoundParticle[] roundParticles = new RoundParticle[2];

    // Start is called before the first frame update
    void Start()
    {
        float u2,tmp;
        plotTexture = new Texture2D(DIM_X,DIM_Y);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),UnityEngine.Vector2.zero);
        
        speed = new float[DIM_X,DIM_Y];
        u = new float[DIM_X,DIM_Y];
        v = new float[DIM_X,DIM_Y];
        fx = new float[DIM_X,DIM_Y];
        fy = new float[DIM_X,DIM_Y];
        rho = new float[DIM_X,DIM_Y];
        f = new float[9,DIM_X,DIM_Y];
        f0 = new float[9,DIM_X,DIM_Y];
        ftmp = new float[9,DIM_X,DIM_Y];

        nu = (tau - 0.5f)/3.0f;

        roundParticles[0] = new RoundParticle(1f,70.0f/200.0f*(float)(DIM_X-1),new float[2]{0.5f*(float)(DIM_X-1),0.5f*(float)(DIM_X-1)});
        roundParticles[1] = new RoundParticle(1f,45.0f/200.0f*(float)(DIM_X-1),new float[2]{0.5f*(float)(DIM_X-1),0.5f*(float)(DIM_X-1)},u0/(45.0f/200.0f*(float)(DIM_X-1)));

        maxSpeed = 0f;
        minSpeed = Mathf.Infinity;
        maxRho = 0f;
        minRho = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u[i,j] = 0f; v[i,j] = 0.0f;
                fx[i,j] = 0.0f; fy[i,j] = 0.0f;
                rho[i,j] = 1.0f;
                u2 = u[i,j]*u[i,j] + v[i,j]*v[i,j];  
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i,j] + cy[k]*v[i,j];  
                    f0[k,i,j] = w[k]*rho[i,j]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
                    f[k,i,j] = f0[k,i,j];
                }
                speed[i,j] = Mathf.Sqrt(u2);
                maxSpeed = Mathf.Max(maxSpeed,speed[i,j]);
                minSpeed = Mathf.Min(minSpeed,speed[i,j]);
                maxRho = Mathf.Max(maxRho,rho[i,j]);
                maxRho = Mathf.Min(maxRho,rho[i,j]);
            } 
        }
    }

    void LBMStep()
    {
        Collision();
        Streaming();
        // BouncebackBoundaries();
        UpdateSpeedAndDensity();
        ImmersedBoundary();
    }

    // Update is called once per frame
    void Update()
    {
        for (int i = 0; i < loopCount; i++)
        {
            LBMStep();
        }
        UpdatePlot();
    }

    void UpdatePlot()
    {
        if(maxSpeed == 0f)maxSpeed = 1f;
        if(maxRho == 0f)maxRho = 1f;
        if(fixHeatMapMax)
        {
            maxSpeed = fixedMaxSpeed;
            maxRho = fixedMaxRho;
        }
        if(fixHeatMapMin)
        {
            minSpeed = fixedMinSpeed;
            minSpeed = fixedMinRho;
        }
        for (int i = 0; i < plotPixels.Length; i++)
        {
            if(mode == HeatMapMode.Speed) 
            plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X,i/DIM_X]-minSpeed,maxSpeed-minSpeed);
            else
            plotPixels[i] = colorHeatMap.GetColorForValue(rho[i%DIM_X,i/DIM_X]-minRho,maxRho-minSpeed);
        }
        for (int n = 0; n < 2; n++)
        {
            // roundParticles[n].PlotParticleFill(ref plotPixels,DIM_X);
            roundParticles[n].PlotParticlePerimeter(ref plotPixels,DIM_X);
        }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
        vectorField.UpdateVectors(u,v,DIM_Y,maxSpeed,vectorFieldOn);
    }

    void Collision()
    {
        float tmp,u2;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u2 = u[i,j]*u[i,j] + v[i,j]*v[i,j];  
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i,j] + cy[k]*v[i,j];  
                    f0[k,i,j] = w[k]*rho[i,j]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
                    f[k,i,j] = f[k,i,j] - (f[k,i,j] - f0[k,i,j])/tau + 3f*w[k]*(fx[i,j]*cx[k] + fy[i,j]*cy[k]);
                }
                // reset force;
                fx[i,j] = 0.0f;
                fy[i,j] = 0.0f;
            }   
        }
    }

    void Streaming()
    {
        ftmp = (float[,,])(f.Clone());

        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            { 
                for(int k = 0; k < 9; k++)
                {
                    // periodic boundary
                    int im = (i + (int)cx[k] + DIM_X)%DIM_X; 
                    int jm = (j + (int)cy[k] + DIM_Y)%DIM_Y;
                    f[k,im,jm] = ftmp[k,i,j];
                    // int im = i + (int)cx[k]; 
                    // int jm = j + (int)cy[k];
                    // if((jm!=DIM_Y&&jm!=-1) && (im!=DIM_X&&im!=-1))
                    // {
                    //     f[k,im,jm] = ftmp[k,i,j];
                    // }
                } 
            }
        }
    }

    void BouncebackBoundaries()
    {
        for (int i = 0; i < DIM_X; i++)
        {
            f[4,i,DIM_Y-1] = f[2,i,DIM_Y-1];
            f[7,i,DIM_Y-1] = f[5,i,DIM_Y-1];
            f[8,i,DIM_Y-1] = f[6,i,DIM_Y-1]; 
            f[2,i,0] = f[4,i,0]; 
            f[5,i,0] = f[7,i,0]; 
            f[6,i,0] = f[8,i,0]; 
        }
        for (int j = 0; j < DIM_Y; j++)
        {
            f[3,DIM_X-1,j] = f[1,DIM_X-1,j];
            f[6,DIM_X-1,j] = f[8,DIM_X-1,j];
            f[7,DIM_X-1,j] = f[5,DIM_X-1,j];
            f[1,0,j] = f[3,0,j]; 
            f[5,0,j] = f[7,0,j]; 
            f[8,0,j] = f[6,0,j]; 
        }
    }

    void UpdateSpeedAndDensity()
    {
        maxSpeed = 0f;
        minSpeed = Mathf.Infinity;
        maxRho = 0f;
        minRho = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                rho[i,j] = f[0,i,j]; 
                u[i,j] = 0; v[i,j] = 0;
                for(int k = 1; k <= 8; k++)
                {
                    rho[i,j] = rho[i,j] + f[k,i,j];
                    u[i,j] =   u[i,j] + f[k,i,j]*cx[k];
                    v[i,j] =   v[i,j] + f[k,i,j]*cy[k];
                } 
                u[i,j] = u[i,j]/rho[i,j];
                v[i,j] = v[i,j]/rho[i,j];
                speed[i,j] = Mathf.Sqrt(u[i,j]*u[i,j] + v[i,j]*v[i,j]);
                maxSpeed = Mathf.Max(maxSpeed,speed[i,j]);
                maxRho = Mathf.Max(maxRho,rho[i,j]);
                minSpeed = Mathf.Min(minSpeed,speed[i,j]);
                minRho = Mathf.Min(minRho,rho[i,j]);
            } 
        }
    }

    void ImmersedBoundary()
    {
        float tmp1,tmp2,tmp3;
        for(int n = 0; n < 2; n++) 
        { 
            for(int m = 0; m < roundParticles[n].perimeterPointCount ; m++) 
            {
                roundParticles[n].perimeterFluidVel[m,0] = 0f;
                roundParticles[n].perimeterFluidVel[m,1] = 0f;
                // 固体表面の速度を計算
                for(int i = (int)roundParticles[n].perimeterPos[m,0] - 3; i < (int)roundParticles[n].perimeterPos[m,0] + 3; i++)
                {
                    for(int j = (int)roundParticles[n].perimeterPos[m,1] - 3; j < (int)roundParticles[n].perimeterPos[m,1] + 3; j++)
                    {
                        tmp1 = Mathf.Abs(roundParticles[n].perimeterPos[m,0] - (float)i);
                        tmp2 = Mathf.Abs(roundParticles[n].perimeterPos[m,1] - (float)j);
                        if(tmp1 <= 2.0f)
                        {
                            tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp1/2.0f))/4.0f;
                        } 
                        else 
                        {
                            tmp3 = 0.0f;
                        }
                        if(tmp2 <= 2.0f)
                        {
                            tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp2/2.0f))/4.0f*tmp3;
                        } 
                        else 
                        {
                            tmp3 = 0.0f;
                        }
                        if((j<DIM_Y&&j>=0) && (i<DIM_X&&i>=0))
                        {
                            roundParticles[n].perimeterFluidVel[m,0] += u[i,j]*tmp3;
                            roundParticles[n].perimeterFluidVel[m,1] += v[i,j]*tmp3;
                        }
                    } 
                }
                roundParticles[n].forceOnPerimeter[m,0] = roundParticles[n].perimeterVel[m,0] - roundParticles[n].perimeterFluidVel[m,0];
                roundParticles[n].forceOnPerimeter[m,1] = roundParticles[n].perimeterVel[m,1] - roundParticles[n].perimeterFluidVel[m,1];
                // 固体が外部に与える力を計算
                for(int i = (int)roundParticles[n].perimeterPos[m,0] - 3; i < (int)roundParticles[n].perimeterPos[m,0] + 3; i++)
                {
                    for(int j = (int)roundParticles[n].perimeterPos[m,1] - 3; j < (int)roundParticles[n].perimeterPos[m,1] + 3; j++)
                    {
                        tmp1 = Mathf.Abs(roundParticles[n].perimeterPos[m,0] - (float)i);
                        tmp2 = Mathf.Abs(roundParticles[n].perimeterPos[m,1] - (float)j);
                        if(tmp1 <= 2.0f)
                        {
                            tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp1/2.0f))/4.0f;
                        } 
                        else 
                        {
                            tmp3 = 0.0f;
                        }
                        if(tmp2 <= 2.0f)
                        {
                            tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp2/2.0f))/4.0f*tmp3;
                        } 
                        else 
                        {
                            tmp3 = 0.0f;
                        }
                        if((j<DIM_Y&&j>=0) && (i<DIM_X&&i>=0))
                        {
                            fx[i,j] += roundParticles[n].forceOnPerimeter[m,0] * tmp3 * 2.0f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;
                            fy[i,j] += roundParticles[n].forceOnPerimeter[m,1] * tmp3 * 2.0f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;
                        }
                    } 
                }
            } 
        }
    }
}
