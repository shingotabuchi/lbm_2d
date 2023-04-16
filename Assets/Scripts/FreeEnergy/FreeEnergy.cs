using System.Collections;
using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class FreeEnergy : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
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
    float[] rhoArray, uvArray, phiArray, cheArray, forceArray;
    float[] fArray;
    float[] gArray;
    
    float tmp,tmp1,tmp2,u2,u3,beta,kap;
    public float tauf =  0.7f;
    public float taug = 0.7f;
    public float gamma = 10.0f, wid = 5.0f, phi0 = 1.0f;
    public float sig = 0.0001f;
    public float radius = 0.25f;

    public int loopCount = 1;

    public float randomBunbo = 10f;

    int DIMSqrd;
    int DIMSqrd9;

    public ComputeShader compute;
    int initMicroVariables,initMacroVariables,collision;
    ComputeBuffer uv,rho,phi;
    ComputeBuffer f,g,che,force;

    public bool nextFrame = false;
    public bool debugMode = true;

    private void OnDestroy() {
        uv.Dispose();
        rho.Dispose();
        phi.Dispose();
        f.Dispose();
        g.Dispose();
        che.Dispose();
        force.Dispose();
    }

    // Start is called before the first frame update
    void Start()
    {
        DIMSqrd = DIM_X * DIM_Y;
        DIMSqrd9 = 9 * DIM_X * DIM_Y;

        plotTexture = new Texture2D(DIM_X,DIM_Y);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),UnityEngine.Vector2.zero);
        ((RectTransform)plotImage.transform).sizeDelta = new Vector2((DIM_X*1080)/DIM_Y,1080);

        uvArray = new float[DIMSqrd*2];
        forceArray = new float[DIMSqrd*2];
        rhoArray = new float[DIMSqrd];
        cheArray = new float[DIMSqrd];
        phiArray = new float[DIMSqrd];
        fArray = new float[DIMSqrd9*2];
        gArray = new float[DIMSqrd9*2];

        beta  = 3.0f/4.0f*sig/wid*Mathf.Pow(phi0,4f);
        kap   = 3.0f/8.0f*sig*wid/Mathf.Pow(phi0,2f);

        uv = new ComputeBuffer(DIMSqrd*2,sizeof(float));
        force = new ComputeBuffer(DIMSqrd*2,sizeof(float));
        rho = new ComputeBuffer(DIMSqrd,sizeof(float));
        phi = new ComputeBuffer(DIMSqrd,sizeof(float));
        che = new ComputeBuffer(DIMSqrd,sizeof(float));
        f = new ComputeBuffer(DIMSqrd9*2,sizeof(float));
        g = new ComputeBuffer(DIMSqrd9*2,sizeof(float));

        maxSpeed = 0f;
        minSpeed = Mathf.Infinity;
        maxPhi = 0f;
        minPhi = Mathf.Infinity;
        maxRho = 0f;
        minRho = Mathf.Infinity;
        maxChe = 0f;
        minChe = Mathf.Infinity;

        compute.SetInt("DIM_X",DIM_X);
        compute.SetInt("DIM_Y",DIM_Y);
        compute.SetInt("DIMSqrd9",DIMSqrd9);
        compute.SetFloat("phi0",phi0);
        compute.SetFloat("wid",wid);
        compute.SetFloat("radius",radius);
        compute.SetFloat("beta",beta);
        compute.SetFloat("kap",kap);
        compute.SetFloat("gamma",gamma);
        compute.SetFloat("tauf",tauf);
        compute.SetFloat("taug",taug);

        initMacroVariables = compute.FindKernel("InitMacroVariables");
        compute.SetBuffer(initMacroVariables,"uv",uv);
        compute.SetBuffer(initMacroVariables,"rho",rho);
        compute.SetBuffer(initMacroVariables,"phi",phi);

        initMicroVariables = compute.FindKernel("InitMicroVariables");
        compute.SetBuffer(initMicroVariables,"uv",uv);
        compute.SetBuffer(initMicroVariables,"rho",rho);
        compute.SetBuffer(initMicroVariables,"phi",phi);
        compute.SetBuffer(initMicroVariables,"f",f);
        compute.SetBuffer(initMicroVariables,"g",g);
        compute.SetBuffer(initMicroVariables,"che",che);

        collision = compute.FindKernel("Collision");
        compute.SetBuffer(collision,"uv",uv);
        compute.SetBuffer(collision,"rho",rho);
        compute.SetBuffer(collision,"phi",phi);
        compute.SetBuffer(collision,"f",f);
        compute.SetBuffer(collision,"g",g);
        compute.SetBuffer(collision,"che",che);
        compute.SetBuffer(collision,"force",force);

        // InitMacro();
        // uv.SetData(uvArray);
        // rho.SetData(rhoArray);
        // phi.SetData(phiArray);
        // f.SetData(fArray);
        // g.SetData(gArray);
        // che.SetData(cheArray);
        compute.Dispatch(initMacroVariables,(DIM_X+7)/8,(DIM_Y+7)/8,1);
        compute.Dispatch(initMicroVariables,(DIM_X+7)/8,(DIM_Y+7)/8,1);
        uv.GetData(uvArray);
        rho.GetData(rhoArray);
        phi.GetData(phiArray);
        f.GetData(fArray);
        g.GetData(gArray);
        che.GetData(cheArray);
        force.GetData(forceArray);
        // Init();
        UpdatePlot();
    }

    
    void InitMacro()
    {
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                uvArray[(i + j*DIM_X)*2 + 0] = 0.0f; uvArray[(i + j*DIM_X)*2 + 1] = 0.0f; 
                rhoArray[i + j*DIM_X] = 1.0f;
                tmp = Mathf.Sqrt(Mathf.Pow(((float)i - (float)(DIM_X-1)*0.5f),2) + Mathf.Pow(((float)j - (float)(DIM_Y-1)*0.5f),2));
                phiArray[i + j*DIM_X] = -phi0*(float)Math.Tanh(2.0f*(tmp - radius*DIM_X)/wid);
            }
        }
    }

    void Init()
    {
        for (int i = 0; i < DIM_X; i++)
        {
            for (int j = 0; j < DIM_Y; j++)
            {
                int ip = (i + 1)%DIM_X;
                int im = (i - 1+DIM_X)%DIM_X;
                int jp = (j + 1)%DIM_Y; 
                int jm = (j - 1+DIM_Y)%DIM_Y; 
                tmp = phiArray[ip + j*DIM_X] + phiArray[im + j*DIM_X] + phiArray[i + jp*DIM_X] + phiArray[i + jm*DIM_X] - 4.0f*phiArray[i + j*DIM_X];
                cheArray[i + j*DIM_X] = 4f*beta*(Mathf.Pow(phiArray[i + j*DIM_X],2) - Mathf.Pow(phi0,2))*phiArray[i + j*DIM_X] - kap*tmp;
                u2 = uvArray[(i + j*DIM_X)*2 + 0] *uvArray[(i + j*DIM_X)*2 + 0] + uvArray[(i + j*DIM_X)*2 + 1] *uvArray[(i + j*DIM_X)*2 + 1] ;
                float f0 = rhoArray[i + j*DIM_X]*(1.0f -3.0f/2.0f*u2 - 15.0f/4.0f*phiArray[i + j*DIM_X]*cheArray[i + j*DIM_X])*4.0f/9.0f;
                float g0 = phiArray[i + j*DIM_X] - 5.0f/9.0f*gamma*cheArray[i + j*DIM_X];

                fArray[(i + j*DIM_X)*9 + 0] = f0;
                gArray[(i + j*DIM_X)*9 + 0] = g0;
                fArray[(i + j*DIM_X)*9 + 0 + DIMSqrd9] = f0;
                gArray[(i + j*DIM_X)*9 + 0 + DIMSqrd9] = g0;

                for (int k = 1; k < 9; k++)
                {
                    tmp = cx[k]*uvArray[(i + j*DIM_X)*2 + 0] + cy[k]*uvArray[(i + j*DIM_X)*2 + 1] ;  
                    f0 = w[k]*rhoArray[i + j*DIM_X]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2 + 3.0f*phiArray[i + j*DIM_X]*cheArray[i + j*DIM_X]);
                    g0 = w[k]*(gamma*cheArray[i + j*DIM_X] + 3.0f*phiArray[i + j*DIM_X]*tmp);
                    fArray[(i + j*DIM_X)*9 + k] = f0;
                    gArray[(i + j*DIM_X)*9 + k] = g0;
                    fArray[(i + j*DIM_X)*9 + k + DIMSqrd9] = f0;
                    gArray[(i + j*DIM_X)*9 + k + DIMSqrd9] = g0;
                }

                maxPhi = Mathf.Max(maxPhi,phiArray[i + j*DIM_X]);
                maxRho = Mathf.Max(maxRho,rhoArray[i + j*DIM_X]);
                maxChe = Mathf.Max(maxChe,cheArray[i + j*DIM_X]);
                minPhi = Mathf.Min(minPhi,phiArray[i + j*DIM_X]);
                minRho = Mathf.Min(minRho,rhoArray[i + j*DIM_X]);
                minChe = Mathf.Min(minChe,cheArray[i + j*DIM_X]);
            }            
        }
    }

    void LBMStep()
    {
        Collision();
        // compute.Dispatch(collision,(DIM_X+7)/8,(DIM_Y+7)/8,1);
        // uv.GetData(uvArray);
        // rho.GetData(rhoArray);
        // phi.GetData(phiArray);
        // f.GetData(fArray);
        // g.GetData(gArray);
        // che.GetData(cheArray);
        // force.GetData(forceArray);

        Streaming();
        // Boundaries();
        UpdateSpeedAndPotentials();

        // uv.SetData(uvArray);
        // rho.SetData(rhoArray);
        // phi.SetData(phiArray);
        // f.SetData(fArray);
        // g.SetData(gArray);
        // che.SetData(cheArray);
        // force.SetData(forceArray);
    }

    // Update is called once per frame
    void Update()
    {
        if(!nextFrame && debugMode) return;
        nextFrame = false;
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
                if(mode == HeatMapMode.ChemicalPotential)
                plotPixels[i] = colorHeatMap.GetColorForValue(cheArray[i%DIM_X+ (i/DIM_X)*DIM_X]-minChe,maxChe-minChe);
                else if(mode == HeatMapMode.OrderParameter)
                plotPixels[i] = colorHeatMap.GetColorForValue(phiArray[i%DIM_X+ (i/DIM_X)*DIM_X]-minPhi,maxPhi - minPhi);
                else
                plotPixels[i] = colorHeatMap.GetColorForValue(rhoArray[i%DIM_X+ (i/DIM_X)*DIM_X]-minRho,maxRho-minRho);
            }
            else
            {
                if(mode == HeatMapMode.ChemicalPotential)
                plotPixels[i] = colorHeatMap.GetColorForValue(cheArray[i%DIM_X+ (i/DIM_X)*DIM_X],maxChe);
                else if(mode == HeatMapMode.OrderParameter)
                plotPixels[i] = colorHeatMap.GetColorForValue(phiArray[i%DIM_X+ (i/DIM_X)*DIM_X],maxPhi);
                else
                plotPixels[i] = colorHeatMap.GetColorForValue(rhoArray[i%DIM_X+ (i/DIM_X)*DIM_X],maxRho);
            }
        }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
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

                u2 = uvArray[(i + j*DIM_X)*2 + 0] *uvArray[(i + j*DIM_X)*2 + 0] + uvArray[(i + j*DIM_X)*2 + 1] *uvArray[(i + j*DIM_X)*2 + 1] ;
                u3 = uvArray[(i + j*DIM_X)*2 + 0] *forceArray[(i + j*DIM_X)*2 + 0]  + uvArray[(i + j*DIM_X)*2 + 1] *forceArray[(i + j*DIM_X)*2 + 1] ;
                forceArray[(i + j*DIM_X)*2 + 0]  = cheArray[i + j*DIM_X] * (phiArray[ip + j*DIM_X] - phiArray[im + j*DIM_X])*0.5f;
                forceArray[(i + j*DIM_X)*2 + 1]  = cheArray[i + j*DIM_X] * (phiArray[i + jp*DIM_X] - phiArray[i + jm*DIM_X])*0.5f;

                float f0 = rhoArray[i + j*DIM_X]*(1.0f -3.0f/2.0f*u2 - 15.0f/4.0f*phiArray[i + j*DIM_X]*cheArray[i + j*DIM_X])*4.0f/9.0f;
                float g0 = phiArray[i + j*DIM_X] - 5.0f/9.0f*gamma*cheArray[i + j*DIM_X];
                float fi = -4.0f/3.0f*(1.0f - 0.5f/tauf)*u3;
                fArray[(i + j*DIM_X)*9 + 0] = fArray[(i + j*DIM_X)*9 + 0 + DIMSqrd9] - (fArray[(i + j*DIM_X)*9 + 0 + DIMSqrd9] - f0)/tauf + fi;
                gArray[(i + j*DIM_X)*9 + 0] = gArray[(i + j*DIM_X)*9 + 0 + DIMSqrd9] - (gArray[(i + j*DIM_X)*9 + 0 + DIMSqrd9] - g0)/taug;

                for (int k = 1; k < 9; k++)
                {
                    tmp = cx[k]*uvArray[(i + j*DIM_X)*2 + 0] + cy[k]*uvArray[(i + j*DIM_X)*2 + 1] ;      
                    tmp1 = cx[k]*forceArray[(i + j*DIM_X)*2 + 0]  + cy[k]*forceArray[(i + j*DIM_X)*2 + 1] ;
                    tmp2 = uvArray[(i + j*DIM_X)*2 + 0] *forceArray[(i + j*DIM_X)*2 + 0] *cx[k]*cx[k]
                            + uvArray[(i + j*DIM_X)*2 + 0] *forceArray[(i + j*DIM_X)*2 + 1] *cx[k]*cy[k]
                            + uvArray[(i + j*DIM_X)*2 + 1] *forceArray[(i + j*DIM_X)*2 + 0] *cy[k]*cx[k]
                            + uvArray[(i + j*DIM_X)*2 + 1] *forceArray[(i + j*DIM_X)*2 + 1] *cy[k]*cy[k];   
                    f0 = w[k]*rhoArray[i + j*DIM_X]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2 + 3.0f*phiArray[i + j*DIM_X]*cheArray[i + j*DIM_X]);
                    g0 = w[k]*(gamma*cheArray[i + j*DIM_X] + 3.0f*phiArray[i + j*DIM_X]*tmp);
                    fi = w[k]*(1.0f - 0.5f/tauf)*(3.0f*tmp1 + 9.0f*tmp2 - 3.0f*u3);
                    fArray[(i + j*DIM_X)*9 + k] = fArray[(i + j*DIM_X)*9 + k + DIMSqrd9] - (fArray[(i + j*DIM_X)*9 + k + DIMSqrd9] - f0)/tauf + fi;
                    gArray[(i + j*DIM_X)*9 + k] = gArray[(i + j*DIM_X)*9 + k + DIMSqrd9] - (gArray[(i + j*DIM_X)*9 + k + DIMSqrd9] - g0)/taug;
                }
            }   
        }
    }

    void Streaming()
    {
        // ftmp = (float[])(fArray.Clone());
        // gtmp = (float[])(gArray.Clone());

        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            { 
                for(int k = 0; k < 9; k++)
                {
                    // periodic boundary
                    int im = (i + (int)cx[k] + DIM_X)%DIM_X; 
                    int jm = (j + (int)cy[k] + DIM_Y)%DIM_Y;
                    fArray[(im + jm*DIM_X)*9 + k + DIMSqrd9] = fArray[(i + j*DIM_X)*9 + k];
                    gArray[(im + jm*DIM_X)*9 + k + DIMSqrd9] = gArray[(i + j*DIM_X)*9 + k];
                    // if((jm!=DIM_Y&&jm!=-1) && (im!=DIM_X&&im!=-1))
                    // {
                    //     fArray[(im + jm*DIM_X)*9 + k] = ftmp[(i + j*DIM_X)*9 + k];
                    //     gArray[(im + jm*DIM_X)*9 + k] = gtmp[(i + j*DIM_X)*9 + k];
                    // }
                } 
            }
        }
    }

    // void Boundaries()
    // {
    //     for (int j = 0; j < DIM_Y; j++)
    //     {
    //         fArray[1,0,j] = fArray[3,0,j];
    //         fArray[5,0,j] = fArray[7,0,j];
    //         fArray[8,0,j] = fArray[6,0,j];
    //         fArray[3,DIM_X-1,j] = fArray[1,DIM_X-1,j]; 
    //         fArray[7,DIM_X-1,j] = fArray[5,DIM_X-1,j]; 
    //         fArray[6,DIM_X-1,j] = fArray[8,DIM_X-1,j]; 

    //         gArray[1,0,j] = gArray[3,0,j];
    //         gArray[5,0,j] = gArray[7,0,j];
    //         gArray[8,0,j] = gArray[6,0,j];
    //         gArray[3,DIM_X-1,j] = gArray[1,DIM_X-1,j];
    //         gArray[7,DIM_X-1,j] = gArray[5,DIM_X-1,j];
    //         gArray[6,DIM_X-1,j] = gArray[8,DIM_X-1,j];
    //     }
    //     for(int i = 0; i < DIM_X; i++){
    //         fArray[4,i,DIM_Y-1] = fArray[2,i,DIM_Y-1];
    //         fArray[7,i,DIM_Y-1] = fArray[5,i,DIM_Y-1];
    //         fArray[8,i,DIM_Y-1] = fArray[6,i,DIM_Y-1]; 
    //         fArray[2,i,0] = fArray[4,i,0]; 
    //         fArray[5,i,0] = fArray[7,i,0]; 
    //         fArray[6,i,0] = fArray[8,i,0]; 
    //         gArray[4,i,DIM_Y-1] = gArray[2,i,DIM_Y-1];
    //         gArray[7,i,DIM_Y-1] = gArray[5,i,DIM_Y-1];
    //         gArray[8,i,DIM_Y-1] = gArray[6,i,DIM_Y-1];
    //         gArray[2,i,0] = gArray[4,i,0]; 
    //         gArray[5,i,0] = gArray[7,i,0]; 
    //         gArray[6,i,0] = gArray[8,i,0]; 
    //     }
    // }

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

                rhoArray[i + j*DIM_X] = fArray[(i + j*DIM_X)*9 + 0 + DIMSqrd9]; 
                uvArray[(i + j*DIM_X)*2 + 0] = 0; uvArray[(i + j*DIM_X)*2 + 1] = 0;
                phiArray[i + j*DIM_X] = gArray[(i + j*DIM_X)*9 + 0 + DIMSqrd9];
                for(int k = 1; k <= 8; k++)
                {
                    rhoArray[i + j*DIM_X] = rhoArray[i + j*DIM_X] + fArray[(i + j*DIM_X)*9 + k + DIMSqrd9];
                    uvArray[(i + j*DIM_X)*2 + 0] =   uvArray[(i + j*DIM_X)*2 + 0] + fArray[(i + j*DIM_X)*9 + k + DIMSqrd9]*cx[k];
                    uvArray[(i + j*DIM_X)*2 + 1] =   uvArray[(i + j*DIM_X)*2 + 1] + fArray[(i + j*DIM_X)*9 + k + DIMSqrd9]*cy[k];
                    phiArray[i + j*DIM_X] =   phiArray[i + j*DIM_X] + gArray[(i + j*DIM_X)*9 + k + DIMSqrd9];
                    
                } 
                uvArray[(i + j*DIM_X)*2 + 0] = uvArray[(i + j*DIM_X)*2 + 0] /rhoArray[i + j*DIM_X];
                uvArray[(i + j*DIM_X)*2 + 1] = uvArray[(i + j*DIM_X)*2 + 1] /rhoArray[i + j*DIM_X];

                uvArray[(i + j*DIM_X)*2 + 0] = uvArray[(i + j*DIM_X)*2 + 0] + forceArray[(i + j*DIM_X)*2 + 0] *0.5f;
                uvArray[(i + j*DIM_X)*2 + 1] = uvArray[(i + j*DIM_X)*2 + 1] + forceArray[(i + j*DIM_X)*2 + 1] *0.5f;
                tmp = phiArray[ip + j*DIM_X] + phiArray[im + j*DIM_X] + phiArray[i + jp*DIM_X] + phiArray[i + jm*DIM_X] - 4.0f*phiArray[i + j*DIM_X];
                cheArray[i + j*DIM_X] = 4f*beta*(Mathf.Pow(phiArray[i + j*DIM_X],2) - Mathf.Pow(phi0,2))*phiArray[i + j*DIM_X] - kap*tmp;

                maxPhi = Mathf.Max(maxPhi,phiArray[i + j*DIM_X]);
                maxRho = Mathf.Max(maxRho,rhoArray[i + j*DIM_X]);
                maxChe = Mathf.Max(maxChe,cheArray[i + j*DIM_X]);
                minPhi = Mathf.Min(minPhi,phiArray[i + j*DIM_X]);
                minRho = Mathf.Min(minRho,rhoArray[i + j*DIM_X]);
                minChe = Mathf.Min(minChe,cheArray[i + j*DIM_X]);
            } 
        }
    }
}
