using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LBMNatConv : MonoBehaviour
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
    float[,,] mf,mf0,mcf;
    float[,,] mg,mg0,mcg;
    float[,] mmu = new float[9,9],ui = new float[9,9],sf = new float[9,9],msf = new float[9,9];
    float[,] mme = new float[5,5],ei = new float[5,5],sg = new float[5,5],msg = new float[5,5];
    
    float umax, umin, tmp, u2, nu, chi, norm, taug, rbetag, h;
    public RTType rtType = RTType.MRT;
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
        mf = new float[9,DIM_X,DIM_Y];
        mf0 = new float[9,DIM_X,DIM_Y];
        mcf = new float[9,DIM_X,DIM_Y];
        mg = new float[5,DIM_X,DIM_Y];
        mg0 = new float[5,DIM_X,DIM_Y];
        mcg = new float[5,DIM_X,DIM_Y];
        InitMatrices();

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

                    if(rtType == RTType.SRT)
                    {
                        f[k,i,j] = f[k,i,j] - (f[k,i,j] - f0[k,i,j])/tauf;
                        if(k<5)
                        g[k,i,j] = g[k,i,j] - (g[k,i,j] - g0[k,i,j])/taug;
                    }
                }
            }   
        }

        if(rtType == RTType.MRT)
        {
            for (int i = 0; i < DIM_X; i++)
            {
                for (int j = 0; j < DIM_Y; j++)
                {
                    for(int k = 0; k <= 8; k++){
                        mf[k,i,j] = 0.0f; mf0[k,i,j] = 0.0f;
                        for(int m = 0; m <= 8; m++){
                            mf[k,i,j] =  mf[k,i,j] + mmu[k,m]* f[m,i,j];
                            mf0[k,i,j] = mf0[k,i,j] + mmu[k,m]*f0[m,i,j];
                        }
                        if(k<5)
                        {
                            mg[k,i,j] = 0.0f; mg0[k,i,j]= 0.0f;
                            for(int m = 0; m <= 4; m++){
                                mg[k,i,j] =  mg[k,i,j] + mme[k,m]* g[m,i,j];
                                mg0[k,i,j] = mg0[k,i,j] + mme[k,m]*g0[m,i,j];
                            }
                        }
                    }
                }
            }
            for (int i = 0; i < DIM_X; i++)
            {
                for (int j = 0; j < DIM_Y; j++)
                {
                    for(int k = 0; k <= 8; k++){
                        mcf[k,i,j] = 0.0f;
                        for(int m = 0; m <= 8; m++){
                            mcf[k,i,j] = mcf[k,i,j] + msf[k,m]*(mf[m,i,j] - mf0[m,i,j]);
                        }
                        if(k<5)
                        {
                            mcg[k,i,j] = 0.0f;
                            for(int m = 0; m <= 4; m++){
                                mcg[k,i,j] = mcg[k,i,j] + msg[k,m]*(mg[m,i,j] - mg0[m,i,j]);
                            }
                        }
                    }
                }
            }
            for (int i = 0; i < DIM_X; i++)
            {
                for (int j = 0; j < DIM_Y; j++)
                {
                    for(int k = 0; k <= 8; k++){
                        f[k,i,j] = f[k,i,j] - mcf[k,i,j];
                        if(k<5)
                        {
                            g[k,i,j] = g[k,i,j] - mcg[k,i,j];
                        }
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

    void InitMatrices()
    {
        //.. collision matrix (f)
        mmu[0,0]= 1.0f; mmu[0,1]= 1.0f; mmu[0,2]= 1.0f;
        mmu[0,3]= 1.0f; mmu[0,4]= 1.0f; mmu[0,5]= 1.0f;
        mmu[0,6]= 1.0f; mmu[0,7]= 1.0f; mmu[0,8]= 1.0f;

        mmu[1,0]=-4.0f; mmu[1,1]=-1.0f; mmu[1,2]=-1.0f;
        mmu[1,3]=-1.0f; mmu[1,4]=-1.0f; mmu[1,5]= 2.0f;
        mmu[1,6]= 2.0f; mmu[1,7]= 2.0f; mmu[1,8]= 2.0f;

        mmu[2,0]= 4.0f; mmu[2,1]=-2.0f; mmu[2,2]=-2.0f;
        mmu[2,3]=-2.0f; mmu[2,4]=-2.0f; mmu[2,5]= 1.0f;
        mmu[2,6]= 1.0f; mmu[2,7]= 1.0f; mmu[2,8]= 1.0f;

        mmu[3,0]= 0.0f; mmu[3,1]= 1.0f; mmu[3,2]= 0.0f;
        mmu[3,3]=-1.0f; mmu[3,4]= 0.0f; mmu[3,5]= 1.0f;
        mmu[3,6]=-1.0f; mmu[3,7]=-1.0f; mmu[3,8]= 1.0f;

        mmu[4,0]= 0.0f; mmu[4,1]=-2.0f; mmu[4,2]= 0.0f;
        mmu[4,3]= 2.0f; mmu[4,4]= 0.0f; mmu[4,5]= 1.0f;
        mmu[4,6]=-1.0f; mmu[4,7]=-1.0f; mmu[4,8]= 1.0f;

        mmu[5,0]= 0.0f; mmu[5,1]= 0.0f; mmu[5,2]= 1.0f;
        mmu[5,3]= 0.0f; mmu[5,4]=-1.0f; mmu[5,5]= 1.0f;
        mmu[5,6]= 1.0f; mmu[5,7]=-1.0f; mmu[5,8]=-1.0f;

        mmu[6,0]= 0.0f; mmu[6,1]= 0.0f; mmu[6,2]=-2.0f;
        mmu[6,3]= 0.0f; mmu[6,4]= 2.0f; mmu[6,5]= 1.0f;
        mmu[6,6]= 1.0f; mmu[6,7]=-1.0f; mmu[6,8]=-1.0f;

        mmu[7,0]= 0.0f; mmu[7,1]= 1.0f; mmu[7,2]=-1.0f;
        mmu[7,3]= 1.0f; mmu[7,4]=-1.0f; mmu[7,5]= 0.0f;
        mmu[7,6]= 0.0f; mmu[7,7]= 0.0f; mmu[7,8]= 0.0f;
        
        mmu[8,0]= 0.0f; mmu[8,1]= 0.0f; mmu[8,2]= 0.0f;
        mmu[8,3]= 0.0f; mmu[8,4]= 0.0f; mmu[8,5]= 1.0f;
        mmu[8,6]=-1.0f; mmu[8,7]= 1.0f; mmu[8,8]=-1.0f;

        //
        ui[0,0]= 1.0f/9.0f ; ui[0,1]= -1.0f/9.0f; ui[0,2]= 1.0f/9.0f ;
        ui[0,3]= 0.0f     ; ui[0,4]= 0.0f     ; ui[0,5]= 0.0f     ;
        ui[0,6]= 0.0f     ; ui[0,7]= 0.0f     ; ui[0,8]= 0.0f     ;

        ui[1,0]= 1.0f/9.0f ; ui[1,1]=-1.0f/36.0f; ui[1,2]=-1.0f/18.0f;
        ui[1,3]= 1.0f/6.0f ; ui[1,4]=-1.0f/6.0f ; ui[1,5]= 0.0f     ;
        ui[1,6]= 0.0f     ; ui[1,7]= 1.0f/4.0f ; ui[1,8]= 0.0f     ;

        ui[2,0]= 1.0f/9.0f ; ui[2,1]=-1.0f/36.0f; ui[2,2]=-1.0f/18.0f;
        ui[2,3]= 0.0f     ; ui[2,4]= 0.0f     ; ui[2,5]= 1.0f/6.0f ;
        ui[2,6]=-1.0f/6.0f ; ui[2,7]=-1.0f/4.0f ; ui[2,8]= 0.0f     ;

        ui[3,0]= 1.0f/9.0f ; ui[3,1]=-1.0f/36.0f; ui[3,2]=-1.0f/18.0f;
        ui[3,3]=-1.0f/6.0f ; ui[3,4]= 1.0f/6.0f ; ui[3,5]= 0.0f     ;
        ui[3,6]= 0.0f     ; ui[3,7]= 1.0f/4.0f ; ui[3,8]= 0.0f     ;

        ui[4,0]= 1.0f/9.0f ; ui[4,1]=-1.0f/36.0f; ui[4,2]=-1.0f/18.0f;
        ui[4,3]= 0.0f     ; ui[4,4]= 0.0f     ; ui[4,5]=-1.0f/6.0f ;
        ui[4,6]= 1.0f/6.0f ; ui[4,7]=-1.0f/4.0f ; ui[4,8]= 0.0f     ;

        ui[5,0]= 1.0f/9.0f ; ui[5,1]= 1.0f/18.0f; ui[5,2]= 1.0f/36.0f;
        ui[5,3]= 1.0f/6.0f ; ui[5,4]= 1.0f/12.0f; ui[5,5]= 1.0f/6.0f ;
        ui[5,6]= 1.0f/12.0f; ui[5,7]= 0.0f     ; ui[5,8]= 1.0f/4.0f ;

        ui[6,0]= 1.0f/9.0f ; ui[6,1]= 1.0f/18.0f; ui[6,2]= 1.0f/36.0f;
        ui[6,3]=-1.0f/6.0f ; ui[6,4]=-1.0f/12.0f; ui[6,5]= 1.0f/6.0f ;
        ui[6,6]= 1.0f/12.0f; ui[6,7]= 0.0f     ; ui[6,8]=-1.0f/4.0f ;

        ui[7,0]= 1.0f/9.0f ; ui[7,1]= 1.0f/18.0f; ui[7,2]= 1.0f/36.0f;
        ui[7,3]=-1.0f/6.0f ; ui[7,4]=-1.0f/12.0f; ui[7,5]=-1.0f/6.0f ;
        ui[7,6]=-1.0f/12.0f; ui[7,7]= 0.0f     ; ui[7,8]= 1.0f/4.0f ;
        
        ui[8,0]= 1.0f/9.0f ; ui[8,1]= 1.0f/18.0f; ui[8,2]= 1.0f/36.0f;
        ui[8,3]= 1.0f/6.0f ; ui[8,4]= 1.0f/12.0f; ui[8,5]=-1.0f/6.0f ;
        ui[8,6]=-1.0f/12.0f; ui[8,7]= 0.0f     ; ui[8,8]=-1.0f/4.0f ;
        
        for(int i = 0; i <= 8; i++){ for(int j = 0; j <= 8; j++){
    sf[i,j] = 0.0f;
  } }

  sf[0,0] = 0.0f; sf[1,1] = 1.5f     ; sf[2,2] = 1.4f     ;
  sf[3,3] = 0.0f; sf[4,4] = 1.5f     ; sf[5,5] = 0.0f     ;
  sf[6,6] = 1.5f; sf[7,7] = 1.0f/tauf; sf[8,8] = 1.0f/tauf;

  for(int i = 0; i <= 8; i++){ for(int j = 0; j <= 8; j++){
    msf[i,j] = 0.0f;
    for(int k = 0; k <= 8; k++){
      msf[i,j] = msf[i,j] + ui[i,k]*sf[k,j];
    }
  } }

//.. collision matrix (g)
  mme[0,0]= 1.0f;mme[0,1]= 1.0f;mme[0,2]= 1.0f;mme[0,3]= 1.0f;mme[0,4]= 1.0f;
  mme[1,0]= 0.0f;mme[1,1]= 1.0f;mme[1,2]= 0.0f;mme[1,3]=-1.0f;mme[1,4]= 0.0f;
  mme[2,0]= 0.0f;mme[2,1]= 0.0f;mme[2,2]= 1.0f;mme[2,3]= 0.0f;mme[2,4]=-1.0f;
  mme[3,0]= 4.0f;mme[3,1]=-1.0f;mme[3,2]=-1.0f;mme[3,3]=-1.0f;mme[3,4]=-1.0f;
  mme[4,0]= 0.0f;mme[4,1]= 1.0f;mme[4,2]=-1.0f;mme[4,3]= 1.0f;mme[4,4]=-1.0f;

  ei[0,0]=0.2f;ei[0,1]= 0.0f;ei[0,2]= 0.0f;ei[0,3]= 0.20f;ei[0,4]= 0.00f;
  ei[1,0]=0.2f;ei[1,1]= 0.5f;ei[1,2]= 0.0f;ei[1,3]=-0.05f;ei[1,4]= 0.25f;
  ei[2,0]=0.2f;ei[2,1]= 0.0f;ei[2,2]= 0.5f;ei[2,3]=-0.05f;ei[2,4]=-0.25f;
  ei[3,0]=0.2f;ei[3,1]=-0.5f;ei[3,2]= 0.0f;ei[3,3]=-0.05f;ei[3,4]= 0.25f;
  ei[4,0]=0.2f;ei[4,1]= 0.0f;ei[4,2]=-0.5f;ei[4,3]=-0.05f;ei[4,4]=-0.25f;

  for(int i = 0; i <= 4; i++){ for(int j = 0; j <= 4; j++){
    sg[i,j]  = 0.0f;
  } }
  sg[0,0] = 0.0f; sg[1,1] = 1.0f/taug; sg[2,2] = 1.0f/taug; sg[3,3] = 1.0f; sg[4,4] = 1.0f;

  for(int i = 0; i <= 4; i++){ for(int j = 0; j <= 4; j++){
     msg[i,j] = 0.0f;
     for(int k = 0; k <= 4; k++){
       msg[i,j] = msg[i,j] + ei[i,k] * sg[k,j];
     }
   } }
    }
}
