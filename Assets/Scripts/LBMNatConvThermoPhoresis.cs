using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LBMNatConvThermoPhoresis : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public VectorField vectorField;
    public HeatMapMode mode = HeatMapMode.Speed;
    public BoundaryType[] wallboundaries = new BoundaryType[]{BoundaryType.Adiabatic,BoundaryType.Adiabatic,BoundaryType.Adiabatic,BoundaryType.Adiabatic};
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
    float[,] rho, u, v, temp, forceFromGravityX, forceFromGravityY,speed, forceFromParticlesX,forceFromParticlesY;
    float[,,] f, f0, ftmp;
    float[,,] g, g0, gtmp;
    float umax, umin, tmp, u2, nu, chi, norm, taug, rbetag, h;
    public float pr = 0.71f;
    public float ra =   10000.0f;
    public float tauf = 0.8f;
    public int loopCount = 1;
    public int particleCount = 25;
    public float particleDensity = 1.25f;
    public float particleRadius = 1.25f;
    public float zeta = 1f;
    public float epsw = 1f;
    public float beta = 1f;
    float[] gravity = new float[2];
    RoundParticle[] roundParticles;
    
    public InitialTemperatureDistribution initTemp = InitialTemperatureDistribution.XGradient;

    public GameObject particlePrefab;
    public Transform particleParent;
    Vector2 plotSizeDelta; 
    float plotScale; 
    [Range(0.0f, 1.0f)]
    public float particleTemp = 1f;
    [Range(0.0f, 1.0f)]
    public float particleLowerTemp = 0f;
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
        forceFromGravityX = new float[DIM_X,DIM_Y];
        forceFromGravityY = new float[DIM_X,DIM_Y];
        forceFromParticlesX = new float[DIM_X,DIM_Y];
        forceFromParticlesY = new float[DIM_X,DIM_Y];
        rho = new float[DIM_X,DIM_Y];
        temp = new float[DIM_X,DIM_Y];
        f = new float[9,DIM_X,DIM_Y];
        f0 = new float[9,DIM_X,DIM_Y];
        ftmp = new float[9,DIM_X,DIM_Y];
        g = new float[5,DIM_X,DIM_Y];
        g0 = new float[5,DIM_X,DIM_Y];
        gtmp = new float[5,DIM_X,DIM_Y];
        gravity[0] = 0f;
        gravity[1] = -rbetag/beta;
        roundParticles = new RoundParticle[particleCount];
        float scale = (2f*particleRadius*plotImage.transform.GetComponent<RectTransform>().sizeDelta.x * plotImage.transform.localScale.x)/(DIM_X*particlePrefab.GetComponent<RectTransform>().sizeDelta.x);
        plotSizeDelta = plotImage.transform.GetComponent<RectTransform>().sizeDelta;
        plotScale = plotImage.transform.localScale.x;
        for (int i = 0; i < particleCount; i++)
        {
            float[] spawnPos =  new float[2]{UnityEngine.Random.Range(particleRadius,DIM_X-particleRadius),UnityEngine.Random.Range(particleRadius,DIM_Y-particleRadius)};
            roundParticles[i] = new RoundParticle(particleDensity,particleRadius,spawnPos,Random.Range(0f,2f*Mathf.PI));
            roundParticles[i].obj = Instantiate(particlePrefab,particleParent);
            roundParticles[i].obj.GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
                                                                                + new Vector2((spawnPos[0]*plotSizeDelta.x)/DIM_X,(spawnPos[1]*plotSizeDelta.y)/DIM_Y))*plotScale;
            roundParticles[i].obj.transform.localScale = new Vector3(scale,scale,1);
        }
        maxTemp = 0f;
        minTemp = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u[i,j] = 0.0f; v[i,j] = 0.0f; 
                rho[i,j] = 1.0f;
                forceFromParticlesX[i,j] = 0f;
                forceFromParticlesY[i,j] = 0f;
                if(initTemp == InitialTemperatureDistribution.AllZero) temp[i,j] = 0f;
                else if(initTemp == InitialTemperatureDistribution.XGradient) temp[i,j] = (float)(DIM_X-1 - i)/(float)(DIM_X-1 - 1);
                else if(initTemp == InitialTemperatureDistribution.YGradient) temp[i,j] = (float)(DIM_Y-1 - j)/(float)(DIM_Y-1 - 1);

                u2 = u[i,j]*u[i,j] + v[i,j]*v[i,j];   
                
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i,j] + cy[k]*v[i,j];     
                    f0[k,i,j] = wf[k]*rho[i,j]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
                    f[k,i,j] = f0[k,i,j];
                    if(k<5)
                    {
                        g0[k,i,j] = wg[k]*temp[i,j]*(1.0f + 3.0f*tmp);
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
        if(maxTemp == 0f)maxTemp = 1f;
        for (int i = 0; i < plotPixels.Length; i++)
        {
            if(normalizeHeatMap)
            {
                if(mode == HeatMapMode.Speed) 
                plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X,i/DIM_X]-minSpeed,maxSpeed-minSpeed);
                else
                plotPixels[i] = colorHeatMap.GetColorForValue(temp[i%DIM_X,i/DIM_X]-minTemp,maxTemp-minTemp);
            }
            else
            {
                if(mode == HeatMapMode.Speed) 
                plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X,i/DIM_X],maxSpeed);
                else
                plotPixels[i] = colorHeatMap.GetColorForValue(temp[i%DIM_X,i/DIM_X],maxTemp);
            }
        }
        for (int i = 0; i < particleCount; i++)
        {
            // roundParticles[i].PlotParticleFill(ref plotPixels,DIM_X);
            roundParticles[i].obj.GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
                                                                                + new Vector2((roundParticles[i].pos[0]*plotSizeDelta.x)/DIM_X,(roundParticles[i].pos[1]*plotSizeDelta.y)/DIM_Y))*plotScale;
            roundParticles[i].obj.transform.rotation = Quaternion.Euler(0,0,((roundParticles[i].theta*180f)/Mathf.PI));
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
                        g0[k,i,j] = wg[k]*temp[i,j]*(1.0f + 3.0f*tmp);
                    }

                    f[k,i,j] = f[k,i,j] - (f[k,i,j] - f0[k,i,j])/tauf;
                    if(k<5)
                    g[k,i,j] = g[k,i,j] - (g[k,i,j] - g0[k,i,j])/taug;
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
                forceFromGravityX[i,j] = 0.0f; forceFromGravityY[i,j] = -beta*gravity[1]*(temp[i,j] - 0.5f);
                for (int k = 0; k < 9; k++)
                {
                    f[k,i,j] = f[k,i,j] + 3f*wf[k]*(cx[k]*( forceFromGravityX[i,j] + forceFromParticlesX[i,j] ) + cy[k]*( forceFromGravityY[i,j] + forceFromParticlesY[i,j] ));
                }
                forceFromParticlesX[i,j] = 0f;
                forceFromParticlesY[i,j] = 0f;
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
        maxTemp = particleTemp;
        // minTemp = particleLowerTemp;
        minTemp = Mathf.Infinity;
        maxRho = 0f;
        minRho = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u[i,j] = 0f; v[i,j] = 0f;
                rho[i,j] = f[0,i,j]; 
                temp[i,j] =  g[0,i,j];
                for(int k = 1; k <= 8; k++)
                {
                    rho[i,j] = rho[i,j] + f[k,i,j];
                    u[i,j] =   u[i,j] + f[k,i,j]*cx[k];
                    v[i,j] =   v[i,j] + f[k,i,j]*cy[k];
                    if(k<5){
                    temp[i,j] =   temp[i,j] + g[k,i,j];}
                    
                } 
                u[i,j] = u[i,j]/rho[i,j];
                v[i,j] = v[i,j]/rho[i,j];
                speed[i,j] = Mathf.Sqrt(u[i,j]*u[i,j] + v[i,j]*v[i,j]);
                
                maxTemp = Mathf.Max(maxTemp,temp[i,j]);
                minTemp = Mathf.Min(minTemp,temp[i,j]);
                maxSpeed = Mathf.Max(maxSpeed,speed[i,j]);
                minSpeed = Mathf.Min(minSpeed,speed[i,j]);
                maxRho = Mathf.Max(maxRho,rho[i,j]);
                minRho = Mathf.Min(minRho,rho[i,j]);
            } 
        }
    }

    void ImmersedBoundary()
    {
        float tmp1,tmp2,tmp3;
        for(int n = 0; n < particleCount; n++) 
        { 
            roundParticles[n].forceFromCollisions[0] = 0f;
            roundParticles[n].forceFromCollisions[1] = 0f;

            for (int k = 0; k < particleCount; k++)
            {
                if(k==n) continue;
                for (int i = 0; i < 2; i++)
                {
                    tmp1 = roundParticles[n].ParticleDistance(roundParticles[k]);
                    if(tmp1 < 2.0f*roundParticles[n].radius + zeta)
                    {
                        roundParticles[n].forceFromCollisions[i] += (roundParticles[n].pos[i] - roundParticles[k].pos[i])*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
                    }
                }
            }

            float wallDistance;
            wallDistance = Mathf.Abs(roundParticles[n].pos[1]);
            if(wallDistance < roundParticles[n].radius) roundParticles[n].pos[1] = roundParticles[n].radius;
            wallDistance = Mathf.Abs(DIM_Y-1-roundParticles[n].pos[1]);
            if(wallDistance < roundParticles[n].radius) roundParticles[n].pos[1] = DIM_Y-1-roundParticles[n].radius;
            wallDistance = Mathf.Abs(roundParticles[n].pos[0]);
            if(wallDistance < roundParticles[n].radius) roundParticles[n].pos[0] =roundParticles[n].radius;
            wallDistance = Mathf.Abs(DIM_X-1-roundParticles[n].pos[0]);
            if(wallDistance < roundParticles[n].radius) roundParticles[n].pos[0] = DIM_X-1-roundParticles[n].radius;

            roundParticles[n].forceFromFluid[0] = 0f;
            roundParticles[n].forceFromFluid[1] = 0f;
            roundParticles[n].torque = 0f;
            for(int m = 0; m < roundParticles[n].perimeterPointCount ; m++) 
            {
                roundParticles[n].perimeterFluidVel[m,0] = 0f;
                roundParticles[n].perimeterFluidVel[m,1] = 0f;
                float angle = Vector2.Angle(
                    new Vector2(roundParticles[n].perimeterPos[m,0] - roundParticles[n].pos[0],roundParticles[n].perimeterPos[m,1] - roundParticles[n].pos[1]),
                    new Vector2(-Mathf.Sin(roundParticles[n].theta),Mathf.Cos(roundParticles[n].theta))
                );
                if((int)roundParticles[n].perimeterPos[m,0] >= 0 && (int)roundParticles[n].perimeterPos[m,0] < DIM_X
                    && (int)roundParticles[n].perimeterPos[m,1] >= 0 && (int)roundParticles[n].perimeterPos[m,1] < DIM_Y
                )
                {
                    if(angle>90f)temp[(int)roundParticles[n].perimeterPos[m,0],(int)roundParticles[n].perimeterPos[m,1]] = particleTemp;
                    // else temp[(int)roundParticles[n].perimeterPos[m,0],(int)roundParticles[n].perimeterPos[m,1]] = particleLowerTemp;
                }
                // else temp[(int)roundParticles[n].perimeterPos[m,0],(int)roundParticles[n].perimeterPos[m,1]] = particleLowerTemp;
                // temp[(int)roundParticles[n].perimeterPos[m,0],(int)roundParticles[n].perimeterPos[m,1]] = 1f;
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
                            forceFromParticlesX[i,j] += roundParticles[n].forceOnPerimeter[m,0] * tmp3 * 2.0f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;
                            forceFromParticlesY[i,j] += roundParticles[n].forceOnPerimeter[m,1] * tmp3 * 2.0f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;
                        }
                    } 
                }
                roundParticles[n].forceFromFluid[0] += roundParticles[n].forceOnPerimeter[m,0];
                roundParticles[n].forceFromFluid[1] += roundParticles[n].forceOnPerimeter[m,1];
                roundParticles[n].torque += roundParticles[n].forceOnPerimeter[m,1] * (roundParticles[n].perimeterPos[m,0] - roundParticles[n].pos[0]) 
                                        - roundParticles[n].forceOnPerimeter[m,0] * (roundParticles[n].perimeterPos[m,1] - roundParticles[n].pos[1]);
            } 

            roundParticles[n].forceFromFluid[0] *= -2f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;  
            roundParticles[n].forceFromFluid[1] *= -2f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;  
            roundParticles[n].torque *= -2f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;  

            roundParticles[n].UpdatePosVel(gravity);
            roundParticles[n].UpdateOmegaTheta();
            roundParticles[n].UpdatePerimeter();
        }
    }
}
