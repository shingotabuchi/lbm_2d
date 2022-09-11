using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LBMNatConvMDArray : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public VectorField vectorField;
    public HeatMapMode mode = HeatMapMode.Speed;
    public BoundaryType[] wallboundaries = new BoundaryType[]{BoundaryType.Bounceback,BoundaryType.Bounceback,BoundaryType.Bounceback,BoundaryType.Bounceback};
    [Range(0.0f, 1.0f)]
    public float wallTemp1 = 1f;
    [Range(0.0f, 1.0f)]
    public float wallTemp2 = 0f;
    [Range(0.0f, 1.0f)]
    public float wallTemp3 = 0f;
    [Range(0.0f, 1.0f)]
    public float wallTemp4 = 0f;
    public bool normalizeHeatMap;
    public int DIM_X = 47;
    public int DIM_Y = 47;

    public float maxTemp,minTemp;
    public float maxSpeed,minSpeed;
    public float maxRho,minRho;

    float[] cx = new float[9]{0, 1,    0,   -1,    0,     1,    -1,    -1,     1};
    float[] cy = new float[9]{0, 0,    1,    0,   -1,     1,     1,    -1,    -1};
    float[] wg = new float[5]{1f/3f,1f/6f,1f/6f,1f/6f,1f/6f};
    float[] wf = new float[9]{4f/9f,1f/9f,1f/9f,1f/9f,1f/9f,1f/36f,1f/36f,1f/36f,1f/36f};
    float[,] rho, u, v, e, fx, fy,speed;
    float[,,] f, f0, ftmp;
    float[,,] g, g0, gtmp;
    
    float umax, umin, tmp, u2, nu, chi, norm, taug, rbetag, h;
    public float pr = 0.71f;
    public float ra =   10000.0f;
    public float tauf = 0.8f;

    public int loopCount = 1;

    // Start is called before the first frame update
    void Start()
    {
        plotTexture = new Texture2D(DIM_X,DIM_Y);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),UnityEngine.Vector2.zero);

        h = (float)(DIM_X-1 - 1);
        nu = (tauf - 0.5f)/3.0f;
        chi = nu/pr;
        taug = 3.0f*chi + 0.5f;

        rbetag = ra*nu*chi/h/h/h;
        speed = new float[DIM_X,DIM_Y];
        u = new float[DIM_X,DIM_Y];
        v = new float[DIM_X,DIM_Y];
        fx = new float[DIM_X,DIM_Y];
        fy = new float[DIM_X,DIM_Y];
        rho = new float[DIM_X,DIM_Y];
        e = new float[DIM_X,DIM_Y];
        f = new float[9,DIM_X,DIM_Y];
        f0 = new float[9,DIM_X,DIM_Y];
        ftmp = new float[9,DIM_X,DIM_Y];
        g = new float[5,DIM_X,DIM_Y];
        g0 = new float[5,DIM_X,DIM_Y];
        gtmp = new float[5,DIM_X,DIM_Y];

        maxTemp = 0f;
        minTemp = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u[i,j] = 0.0f; v[i,j] = 0.0f; 
                rho[i,j] = 1.0f;
                e[i,j] = (float)(DIM_X-1 - i)/(float)(DIM_X-1 - 1);

                u2 = u[i,j]*u[i,j] + v[i,j]*v[i,j];   
                
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i,j] + cy[k]*v[i,j];     
                    f0[k,i,j] = wf[k]*rho[i,j]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
                    f[k,i,j] = f0[k,i,j];
                    if(k<5)
                    {
                        g0[k,i,j] = wg[k]*e[i,j]*(1.0f + 3.0f*tmp);
                        g[k,i,j] = g0[k,i,j];
                    }
                }
            }
        }
    }

    void LBMStep()
    {
        Collision();
        Force();
        Streaming();
        Boundaries();
        UpdateSpeedAndTemperature();
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
                else
                plotPixels[i] = colorHeatMap.GetColorForValue(e[i%DIM_X,i/DIM_X]-minTemp,maxTemp-minTemp);
            }
            else
            {
                if(mode == HeatMapMode.Speed) 
                plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X,i/DIM_X],maxSpeed);
                else
                plotPixels[i] = colorHeatMap.GetColorForValue(e[i%DIM_X,i/DIM_X],maxTemp);
            }
        }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
        vectorField.UpdateVectors(u,v,DIM_Y,maxSpeed);
    }

    void Collision()
    {
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u2 = u[i,j]*u[i,j] + v[i,j]*v[i,j];   
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i,j] + cy[k]*v[i,j];     
                    f0[k,i,j] = wf[k]*rho[i,j]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
                    if(k<5)
                    {
                        g0[k,i,j] = wg[k]*e[i,j]*(1.0f + 3.0f*tmp);
                    }

                    f[k,i,j] = f[k,i,j] - (f[k,i,j] - f0[k,i,j])/tauf;
                    if(k<5)
                    {
                        g[k,i,j] = g[k,i,j] - (g[k,i,j] - g0[k,i,j])/taug;
                        gtmp[k,i,j] = g[k,i,j];
                    }
                }
            }   
        }
    }

    void Force()
    {
        for (int i = 0; i < DIM_X; i++)
        {
            for (int j = 0; j < DIM_Y; j++)
            {
                fx[i,j] = 0.0f; fy[i,j] = rbetag*(e[i,j] - 0.5f);
                for (int k = 0; k < 9; k++)
                {
                    f[k,i,j] = f[k,i,j] + 3f*wf[k]*(cx[k]*fx[i,j] + cy[k]*fy[i,j]);
                    ftmp[k,i,j] = f[k,i,j];
                }
            }
        }
    }
    
    void Streaming()
    {
        // ftmp = (float[,,])(f.Clone());
        // gtmp = (float[,,])(g.Clone());

        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            { 
                for(int k = 0; k < 9; k++)
                {
                    int im = i + (int)cx[k]; 
                    int jm = j + (int)cy[k];
                    if((jm!=DIM_Y&&jm!=-1) && (im!=DIM_X&&im!=-1))
                    {
                        f[k,im,jm] = ftmp[k,i,j];
                        if(k<5)
                        g[k,im,jm] = gtmp[k,i,j];
                    }
                } 
            }
        }
    }

    void Boundaries()
    {
        for (int j = 0; j < DIM_Y; j++)
        // for (int j = 1; j < DIM_Y-1; j++)
        {
            f[1,0,j] = f[3,0,j];
            f[5,0,j] = f[7,0,j];
            f[8,0,j] = f[6,0,j];
            f[3,DIM_X-1,j] = f[1,DIM_X-1,j]; 
            f[7,DIM_X-1,j] = f[5,DIM_X-1,j]; 
            f[6,DIM_X-1,j] = f[8,DIM_X-1,j]; 
            
            switch (wallboundaries[0])
            {
                case BoundaryType.Constant: 
                g[0,0,j] = wg[0]*wallTemp1;
                g[1,0,j] = 2f*wg[1]*wallTemp1 - g[2,0,j];
                g[3,0,j] = 2f*wg[1]*wallTemp1 - g[4,0,j];
                break;
                case BoundaryType.Adiabatic:
                for (int k = 0; k < 5; k++)
                {
                    g[k,0,j] = g[k,1,j];
                }
                break;
                case BoundaryType.Bounceback:
                g[1,0,j] = g[3, 0,j ];
                break;
            }
            switch (wallboundaries[1])
            {
                case BoundaryType.Constant: 
                g[0,DIM_X-1,j] = wg[0]*wallTemp2;
                g[2,DIM_X-1,j] = 2f*wg[1]*wallTemp2 - g[1,DIM_X-1,j];
                g[3,DIM_X-1,j] = 2f*wg[1]*wallTemp2 - g[4,DIM_X-1,j];
                break;
                case BoundaryType.Adiabatic:
                for (int k = 0; k < 5; k++)
                {
                    g[k,DIM_X-1,j] = g[k,DIM_X-2,j];
                }
                break;
                case BoundaryType.Bounceback:
                g[3,DIM_X-1,j] =g[1,DIM_X-1,j ];
                break;
            }
            // int jm = j - 1; int jp = j + 1;
            // f[1,   1,j] = f[3,0,j ];
            // f[5,   1,j] = f[7,0,jm];
            // f[8,   1,j] = f[6,0,jp];

            // f[3,DIM_X-1-1,j] = f[1,DIM_X-1,j ];
            // f[7,DIM_X-1-1,j] = f[5,DIM_X-1,jp];
            // f[6,DIM_X-1-1,j] = f[8,DIM_X-1,jm];

            // g[1,   1,j] =-g[3, 0,j ] + 1.0f/3.0f;
            // g[3,DIM_X-1-1,j] =-g[1,DIM_X-1,j ];
        }
        for(int i = 0; i < DIM_X; i++){
        // for(int i = 1; i < DIM_X-1; i++){
            f[4,i,DIM_Y-1] = f[2,i,DIM_Y-1];
            f[7,i,DIM_Y-1] = f[5,i,DIM_Y-1];
            f[8,i,DIM_Y-1] = f[6,i,DIM_Y-1]; 
            f[2,i,0] = f[4,i,0]; 
            f[5,i,0] = f[7,i,0]; 
            f[6,i,0] = f[8,i,0]; 

            // int im = i - 1; int ip = i + 1;
            // f[2,i,   1] = f[4,i , 0];
            // f[5,i,   1] = f[7,im, 0];
            // f[6,i,   1] = f[8,ip, 0];

            // f[4,i,DIM_Y-1-1] = f[2,i ,DIM_Y-1];
            // f[7,i,DIM_Y-1-1] = f[5,ip,DIM_Y-1];
            // f[8,i,DIM_Y-1-1] = f[6,im,DIM_Y-1];

            // g[2,i,   1] = g[4,i , 0];
            // g[4,i,DIM_Y-1-1] = g[2,i ,DIM_Y-1];

            switch (wallboundaries[2])
            {
                case BoundaryType.Constant: 
                g[0,i,DIM_Y-1] = wg[0]*wallTemp3;
                g[1,i,DIM_Y-1] = 2f*wg[1]*wallTemp3 - g[2,i,DIM_Y-1];
                g[4,i,DIM_Y-1] = 2f*wg[1]*wallTemp3 - g[3,i,DIM_Y-1];
                break;
                case BoundaryType.Adiabatic:
                for (int k = 0; k < 5; k++)
                {
                    g[k,i,DIM_Y-1] = g[k,i,DIM_Y-2];
                }
                break;
                case BoundaryType.Bounceback:
                g[4,i,DIM_Y-1] = g[2,i ,DIM_Y-1];
                break;
            }
            switch (wallboundaries[3])
            {
                case BoundaryType.Constant: 
                g[0,i,0] = wg[0]*wallTemp4;
                g[1,i,0] = 2f*wg[1]*wallTemp4 - g[2,i,0];
                g[3,i,0] = 2f*wg[1]*wallTemp4 - g[4,i,0];
                break;
                case BoundaryType.Adiabatic:
                for (int k = 0; k < 5; k++)
                {
                    g[k,i,0] = g[k,i,1];
                }
                break;
                case BoundaryType.Bounceback:
                g[2,i,0] = g[4,i , 0];
                break;
            }
        }
    }

    void UpdateSpeedAndTemperature()
    {
        maxSpeed = 0f;
        minSpeed = Mathf.Infinity;
        maxTemp = 0f;
        minTemp = Mathf.Infinity;
        maxRho = 0f;
        minRho = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u[i,j] = 0f; v[i,j] = 0f;
                rho[i,j] = f[0,i,j]; 
                e[i,j] =  g[0,i,j];
                for(int k = 1; k <= 8; k++)
                {
                    rho[i,j] = rho[i,j] + f[k,i,j];
                    u[i,j] =   u[i,j] + f[k,i,j]*cx[k];
                    v[i,j] =   v[i,j] + f[k,i,j]*cy[k];
                    if(k<5){
                    e[i,j] =   e[i,j] + g[k,i,j];}
                    
                } 
                u[i,j] = u[i,j]/rho[i,j];
                v[i,j] = v[i,j]/rho[i,j];
                speed[i,j] = Mathf.Sqrt(u[i,j]*u[i,j] + v[i,j]*v[i,j]);
                
                maxTemp = Mathf.Max(maxTemp,e[i,j]);
                minTemp = Mathf.Min(minTemp,e[i,j]);
                maxSpeed = Mathf.Max(maxSpeed,speed[i,j]);
                minSpeed = Mathf.Min(minSpeed,speed[i,j]);
                maxRho = Mathf.Max(maxRho,rho[i,j]);
                minRho = Mathf.Min(minRho,rho[i,j]);
            } 
        }
    }
}