using System.Collections;
using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LBMLap : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public VectorField vectorField;
    public bool vectorFieldOn;
    public HeatMapMode mode = HeatMapMode.Speed;
    public bool normalizeHeatMap;
    public int DIM_X = 51;
    public int DIM_Y = 51;

    public float maxRho,minRho;
    public float maxSpeed,minSpeed;
    public float maxPhi,minPhi;
    public float maxChe,minChe;

    float[] cx = new float[9]{0, 1,    0,   -1,    0,     1,    -1,    -1,     1};
    float[] cy = new float[9]{0, 0,    1,    0,   -1,     1,     1,    -1,    -1};
    float[] w = new float[9]{4f/9f,1f/9f,1f/9f,1f/9f,1f/9f,1f/36f,1f/36f,1f/36f,1f/36f};
    float[,] rho, u, v, speed, phi, che, fx ,fy;
    float[,,] f, f0, ftmp, fi;
    float[,,] g, g0, gtmp;
    
    float tmp,tmp1,tmp2,u2,u3,beta,kap;
    public float tauf =  0.7f;
    public float taug = 0.7f;
    public float gamma = 10.0f, wid = 5.0f, phi0 = 1.0f;
    public float sig = 0.0001f;
    public float radius = 12f;

    public int loopCount = 1;

    // Start is called before the first frame update
    void Start()
    {
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
        che = new float[DIM_X,DIM_Y];
        phi = new float[DIM_X,DIM_Y];
        f = new float[9,DIM_X,DIM_Y];
        f0 = new float[9,DIM_X,DIM_Y];
        ftmp = new float[9,DIM_X,DIM_Y];
        fi = new float[9,DIM_X,DIM_Y];
        g = new float[9,DIM_X,DIM_Y];
        g0 = new float[9,DIM_X,DIM_Y];
        gtmp = new float[9,DIM_X,DIM_Y];

        beta  = 3.0f/4.0f*sig/wid*Mathf.Pow(phi0,4f);
        kap   = 3.0f/8.0f*sig*wid/Mathf.Pow(phi0,2f);

        maxSpeed = 0f;
        minSpeed = Mathf.Infinity;
        maxPhi = 0f;
        minPhi = Mathf.Infinity;
        maxRho = 0f;
        minRho = Mathf.Infinity;
        maxChe = 0f;
        minChe = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u[i,j] = 0.0f; v[i,j] = 0.0f; 
                rho[i,j] = 1.0f;
                tmp = Mathf.Sqrt(Mathf.Pow(((float)i - (float)(DIM_X-1)*0.5f),2) + Mathf.Pow(((float)j - (float)(DIM_Y-1)*0.5f),2));
                phi[i,j] = -phi0*(float)Math.Tanh(2.0f*(tmp - radius)/wid);
            }
        }

        for (int i = 0; i < DIM_X; i++)
        {
            for (int j = 0; j < DIM_Y; j++)
            {
                int ip = (i + 1)%DIM_X;
                int im = (i - 1+DIM_X)%DIM_X;
                int jp = (j + 1)%DIM_Y; 
                int jm = (j - 1+DIM_Y)%DIM_Y; 
                tmp = phi[ip,j] + phi[im,j] + phi[i,jp] + phi[i,jm] - 4.0f*phi[i,j];
                che[i,j] = 4f*beta*(Mathf.Pow(phi[i,j],2) - Mathf.Pow(phi0,2))*phi[i,j] - kap*tmp;
                u2 = u[i,j]*u[i,j] + v[i,j]*v[i,j];
                speed[i,j] = Mathf.Sqrt(u2);
                f0[0,i,j] = rho[i,j]*(1.0f -3.0f/2.0f*u2 - 15.0f/4.0f*phi[i,j]*che[i,j])*4.0f/9.0f;
                g0[0,i,j] = phi[i,j] - 5.0f/9.0f*gamma*che[i,j];
                f[0,i,j] = f0[0,i,j];
                g[0,i,j] = g0[0,i,j];
                for (int k = 1; k < 9; k++)
                {
                    tmp = cx[k]*u[i,j] + cy[k]*v[i,j];  
                    f0[k,i,j] = w[k]*rho[i,j]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2 + 3.0f*phi[i,j]*che[i,j]);
                    g0[k,i,j] = w[k]*(gamma*che[i,j] + 3.0f*phi[i,j]*tmp);
                    f[k,i,j] = f0[k,i,j];
                    g[k,i,j] = g0[k,i,j];
                }

                maxSpeed = Mathf.Max(maxSpeed,speed[i,j]);
                maxPhi = Mathf.Max(maxPhi,phi[i,j]);
                maxRho = Mathf.Max(maxRho,rho[i,j]);
                maxChe = Mathf.Max(maxChe,che[i,j]);
                minSpeed = Mathf.Min(minSpeed,speed[i,j]);
                minPhi = Mathf.Min(minPhi,phi[i,j]);
                minRho = Mathf.Min(minRho,rho[i,j]);
                minChe = Mathf.Min(minChe,che[i,j]);
            }            
        }
    }

    void LBMStep()
    {
        Collision();
        Streaming();
        Boundaries();
        UpdateSpeedAndPotentials();
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
        for (int i = 0; i < plotPixels.Length; i++)
        {
            if(normalizeHeatMap)
            {
                if(mode == HeatMapMode.Speed) 
                plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X,i/DIM_X]-minSpeed,maxSpeed-minSpeed);
                else if(mode == HeatMapMode.ChemicalPotential)
                plotPixels[i] = colorHeatMap.GetColorForValue(che[i%DIM_X,i/DIM_X]-minChe,maxChe-minChe);
                else if(mode == HeatMapMode.OrderParameter)
                plotPixels[i] = colorHeatMap.GetColorForValue(phi[i%DIM_X,i/DIM_X]-minPhi,maxPhi-minPhi);
                else
                plotPixels[i] = colorHeatMap.GetColorForValue(rho[i%DIM_X,i/DIM_X]-minRho,maxRho-minRho);
            }
            else
            {
                if(mode == HeatMapMode.Speed) 
                plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X,i/DIM_X],maxSpeed);
                else if(mode == HeatMapMode.ChemicalPotential)
                plotPixels[i] = colorHeatMap.GetColorForValue(che[i%DIM_X,i/DIM_X],maxChe);
                else if(mode == HeatMapMode.OrderParameter)
                plotPixels[i] = colorHeatMap.GetColorForValue(phi[i%DIM_X,i/DIM_X],maxPhi);
                else
                plotPixels[i] = colorHeatMap.GetColorForValue(rho[i%DIM_X,i/DIM_X],maxRho);
            }
        }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
        vectorField.UpdateVectors(u,v,DIM_Y,maxSpeed,vectorFieldOn);
    }

    void Collision()
    {
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                int ip = (i + 1)%DIM_X;
                int im = (i - 1+DIM_X)%DIM_X;
                int jp = (j + 1)%DIM_Y; 
                int jm = (j - 1+DIM_Y)%DIM_Y;

                u2 = u[i,j]*u[i,j] + v[i,j]*v[i,j];
                u3 = u[i,j]*fx[i,j] + v[i,j]*fy[i,j];
                fx[i,j] = che[i,j] * (phi[ip,j ] - phi[im,j ])*0.5f;
                fy[i,j] = che[i,j] * (phi[i ,jp] - phi[i ,jm])*0.5f;

                f0[0,i,j] = rho[i,j]*(1.0f -3.0f/2.0f*u2 - 15.0f/4.0f*phi[i,j]*che[i,j])*4.0f/9.0f;
                g0[0,i,j] = phi[i,j] - 5.0f/9.0f*gamma*che[i,j];
                fi[0,i,j] = -4.0f/3.0f*(1.0f - 0.5f/tauf)*u3;
                f[0,i,j] = f[0,i,j] - (f[0,i,j] - f0[0,i,j])/tauf + fi[0,i,j];
                g[0,i,j] = g[0,i,j] - (g[0,i,j] - g0[0,i,j])/taug;
                for (int k = 1; k < 9; k++)
                {
                    tmp = cx[k]*u[i,j] + cy[k]*v[i,j];      
                    tmp1 = cx[k]*fx[i,j] + cy[k]*fy[i,j];
                    tmp2 = u[i,j]*fx[i,j]*cx[k]*cx[k]
                            + u[i,j]*fy[i,j]*cx[k]*cy[k]
                            + v[i,j]*fx[i,j]*cy[k]*cx[k]
                            + v[i,j]*fy[i,j]*cy[k]*cy[k];   
                    f0[k,i,j] = w[k]*rho[i,j]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2 + 3.0f*phi[i,j]*che[i,j]);
                    g0[k,i,j] = w[k]*(gamma*che[i,j] + 3.0f*phi[i,j]*tmp);
                    fi[k,i,j] = w[k]*(1.0f - 0.5f/tauf)*(3.0f*tmp1 + 9.0f*tmp2 - 3.0f*u3);
                    f[k,i,j] = f[k,i,j] - (f[k,i,j] - f0[k,i,j])/tauf + fi[k,i,j];
                    g[k,i,j] = g[k,i,j] - (g[k,i,j] - g0[k,i,j])/taug;
                }
            }   
        }
    }

    void Streaming()
    {
        ftmp = (float[,,])(f.Clone());
        gtmp = (float[,,])(g.Clone());

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
                    g[k,im,jm] = gtmp[k,i,j];
                    // if((jm!=DIM_Y&&jm!=-1) && (im!=DIM_X&&im!=-1))
                    // {
                    //     f[k,im,jm] = ftmp[k,i,j];
                    //     g[k,im,jm] = gtmp[k,i,j];
                    // }
                } 
            }
        }
    }

    void Boundaries()
    {
        for (int j = 0; j < DIM_Y; j++)
        {
            f[1,0,j] = f[3,0,j];
            f[5,0,j] = f[7,0,j];
            f[8,0,j] = f[6,0,j];
            f[3,DIM_X-1,j] = f[1,DIM_X-1,j]; 
            f[7,DIM_X-1,j] = f[5,DIM_X-1,j]; 
            f[6,DIM_X-1,j] = f[8,DIM_X-1,j]; 

            g[1,0,j] = g[3,0,j];
            g[5,0,j] = g[7,0,j];
            g[8,0,j] = g[6,0,j];
            g[3,DIM_X-1,j] = g[1,DIM_X-1,j];
            g[7,DIM_X-1,j] = g[5,DIM_X-1,j];
            g[6,DIM_X-1,j] = g[8,DIM_X-1,j];
        }
        for(int i = 0; i < DIM_X; i++){
            f[4,i,DIM_Y-1] = f[2,i,DIM_Y-1];
            f[7,i,DIM_Y-1] = f[5,i,DIM_Y-1];
            f[8,i,DIM_Y-1] = f[6,i,DIM_Y-1]; 
            f[2,i,0] = f[4,i,0]; 
            f[5,i,0] = f[7,i,0]; 
            f[6,i,0] = f[8,i,0]; 
            g[4,i,DIM_Y-1] = g[2,i,DIM_Y-1];
            g[7,i,DIM_Y-1] = g[5,i,DIM_Y-1];
            g[8,i,DIM_Y-1] = g[6,i,DIM_Y-1];
            g[2,i,0] = g[4,i,0]; 
            g[5,i,0] = g[7,i,0]; 
            g[6,i,0] = g[8,i,0]; 
        }
    }

    void UpdateSpeedAndPotentials()
    {
        maxSpeed = 0f;
        minSpeed = Mathf.Infinity;
        maxPhi = 0f;
        minPhi = Mathf.Infinity;
        maxRho = 0f;
        minRho = Mathf.Infinity;
        maxChe = 0f;
        minChe = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                int ip = (i + 1)%DIM_X;
                int im = (i - 1+DIM_X)%DIM_X;
                int jp = (j + 1)%DIM_Y; 
                int jm = (j - 1+DIM_Y)%DIM_Y;

                rho[i,j] = f[0,i,j]; 
                u[i,j] = 0; v[i,j] = 0;
                phi[i,j] = g[0,i,j];
                for(int k = 1; k <= 8; k++)
                {
                    rho[i,j] = rho[i,j] + f[k,i,j];
                    u[i,j] =   u[i,j] + f[k,i,j]*cx[k];
                    v[i,j] =   v[i,j] + f[k,i,j]*cy[k];
                    phi[i,j] =   phi[i,j] + g[k,i,j];
                    
                } 
                u[i,j] = u[i,j]/rho[i,j];
                v[i,j] = v[i,j]/rho[i,j];

                u[i,j] = u[i,j] + fx[i,j]*0.5f;
                v[i,j] = v[i,j] + fy[i,j]*0.5f;
                tmp = phi[ip,j] + phi[im,j] + phi[i,jp] + phi[i,jm] - 4.0f*phi[i,j];
                che[i,j] = 4f*beta*(Mathf.Pow(phi[i,j],2) - Mathf.Pow(phi0,2))*phi[i,j] - kap*tmp;

                speed[i,j] = Mathf.Sqrt(u[i,j]*u[i,j] + v[i,j]*v[i,j]);
                maxSpeed = Mathf.Max(maxSpeed,speed[i,j]);
                maxPhi = Mathf.Max(maxPhi,phi[i,j]);
                maxRho = Mathf.Max(maxRho,rho[i,j]);
                maxChe = Mathf.Max(maxChe,che[i,j]);
                minSpeed = Mathf.Min(minSpeed,speed[i,j]);
                minPhi = Mathf.Min(minPhi,phi[i,j]);
                minRho = Mathf.Min(minRho,rho[i,j]);
                minChe = Mathf.Min(minChe,che[i,j]);
            } 
        }
    }
}
