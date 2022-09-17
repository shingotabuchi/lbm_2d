using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LBMVicsek : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public VectorField vectorField;
    public HeatMapMode mode = HeatMapMode.Speed;
    public bool normalizeHeatMap;
    public float tau = 0.56f;
    float tmp,u2;
    float[,] u,v,rho,speed,forceFromParticlesX,forceFromParticlesY;
    float[,,] f,f0,ftmp;
    public float maxSpeed,minSpeed;
    public float maxRho,minRho;
    float[] cx = new float[9]{0,        1,    0,   -1,    0,     1,    -1,    -1,     1};
    float[] cy = new float[9]{0,        0,    1,    0,   -1,     1,     1,    -1,    -1};
    float[] w = new  float[9]{4f/9f,1f/9f,1f/9f,1f/9f,1f/9f,1f/36f,1f/36f,1f/36f,1f/36f};
    public int loopCount = 1;
    
    public GameObject particlePrefab;
    public Transform particleParent;
    public int particleCount;
    RoundParticle[] roundParticles;
    public float particleDensity = 1.25f;
    public float particleRadius = 1.25f;
    public float zeta = 0.1f;
    public float epsw = 10f;
    public int DIM_X,DIM_Y;
    public float particlePropulsion = 0.0001f;
    public float particleViewRange = 2f;
    public float particleFov = 120f;
    float plotScale;
    public float noiseAngleDeg = 30f;
    public float gravity = 0f;

    public bool align;
    Vector2 plotSizeDelta;
    
    // Start is called before the first frame update
    void Start()
    {
        plotTexture = new Texture2D(DIM_X,DIM_Y);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),Vector2.zero);

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

        u = new float[DIM_X,DIM_Y];
        v = new float[DIM_X,DIM_Y];
        rho = new float[DIM_X,DIM_Y];
        speed = new float[DIM_X,DIM_Y];
        f = new float[9,DIM_X,DIM_Y];
        f0 = new float[9,DIM_X,DIM_Y];
        ftmp = new float[9,DIM_X,DIM_Y];
        forceFromParticlesX = new float[DIM_X,DIM_Y];
        forceFromParticlesY = new float[DIM_X,DIM_Y];

        maxSpeed = 0f;
        minSpeed = Mathf.Infinity;
        maxRho = 0f;
        minRho = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u[i,j] = 0.0f;
                v[i,j] = 0.0f;
                u2 = u[i,j]*u[i,j] + v[i,j]*v[i,j];
                rho[i,j] = 1.0f;
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
                minRho = Mathf.Min(minRho,rho[i,j]);
            }
        }
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
        for (int i = 0; i < plotPixels.Length; i++)
        {
            if(normalizeHeatMap)
            {
                if(mode == HeatMapMode.Speed) 
                plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X,i/DIM_X]-minSpeed,maxSpeed-minSpeed);
                else
                plotPixels[i] = colorHeatMap.GetColorForValue(rho[i%DIM_X,i/DIM_X]-minRho,maxRho-minRho);
            }
            else
            {
                if(mode == HeatMapMode.Speed) 
                plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X,i/DIM_X],maxSpeed);
                else
                plotPixels[i] = colorHeatMap.GetColorForValue(rho[i%DIM_X,i/DIM_X],maxRho);
            }
        }
        for (int i = 0; i < particleCount; i++)
        {
            roundParticles[i].obj.GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
                                                                                + new Vector2((roundParticles[i].pos[0]*plotSizeDelta.x)/DIM_X,(roundParticles[i].pos[1]*plotSizeDelta.y)/DIM_Y))*plotScale;
            roundParticles[i].obj.transform.rotation = Quaternion.Euler(0,0,(roundParticles[i].theta*180f)/Mathf.PI);
        }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
        vectorField.UpdateVectors(u,v,DIM_Y,maxSpeed);
    }

    void LBMStep()
    {
        Collision();
        Streaming();
        Boundaries();
        UpdateSpeedAndDensity();
        ImmersedBoundary();
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
                    f0[k,i,j] = w[k]*rho[i,j]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
                    f[k,i,j] += -(f[k,i,j] - f0[k,i,j])/tau;
                    // Force
                    f[k,i,j] += 3f*w[k]*rho[i,j]*(cx[k]*forceFromParticlesX[i,j] + cy[k]*forceFromParticlesY[i,j]);
                    // f[k,i,j] += 3f*w[k]*rho[i,j]*(cx[k]*0.0001f + cy[k]*forceFromParticlesY[i,j]);
                    for (int n = 0; n < particleCount; n++)
                    {
                        if(roundParticles[n].PointIsInParticle(new int[2]{i,j}))
                        {
                            Vector2 force = new Vector2(-Mathf.Sin(roundParticles[n].theta),Mathf.Cos(roundParticles[n].theta));
                            force *= particlePropulsion;
                            f[k,i,j] += 3f*w[k]*rho[i,j]*(cx[k]*(force[0]) + cy[k]*(force[1]));
                            break;
                        }
                    }
                    // if(forceFromParticlesX[i,j]!=0f) print(f[k,i,j]);
                    forceFromParticlesX[i,j] = 0f;
                    forceFromParticlesY[i,j] = 0f;
                } 
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
                for(int k = 0; k <= 8; k++)
                {
                    int im = i + (int)cx[k]; 
                    int jm = j + (int)cy[k];
                    if((jm!=DIM_Y&&jm!=-1) && (im!=DIM_X&&im!=-1))
                    {
                        f[k,im,jm] = ftmp[k,i,j];
                    }
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
        }
        for(int i = 0; i < DIM_X; i++)
        {
            f[4,i,DIM_Y-1] = f[2,i,DIM_Y-1];
            f[7,i,DIM_Y-1] = f[5,i,DIM_Y-1];
            f[8,i,DIM_Y-1] = f[6,i,DIM_Y-1]; 
            f[2,i,0] = f[4,i,0]; 
            f[5,i,0] = f[7,i,0]; 
            f[6,i,0] = f[8,i,0]; 
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
                u[i,j] = 0f; v[i,j] = 0f;
                rho[i,j] = f[0,i,j]; 
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
            // tmp1 = Mathf.Abs(roundParticles[n].pos[1] + roundParticles[n].radius); 
            // if(tmp1 < 2.0f*roundParticles[n].radius + zeta){
            //     roundParticles[n].forceFromCollisions[1] = (roundParticles[n].pos[1] + roundParticles[n].radius)*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
            // }
            // tmp1 = Mathf.Abs(DIM_Y-1-roundParticles[n].pos[1] + roundParticles[n].radius); 
            // if(tmp1 < 2.0f*roundParticles[n].radius + zeta){
            //     roundParticles[n].forceFromCollisions[1] = -(DIM_Y-1-roundParticles[n].pos[1] + roundParticles[n].radius)*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
            // }
            // tmp1 = Mathf.Abs(roundParticles[n].pos[0] + roundParticles[n].radius); 
            // if(tmp1 < 2.0f*roundParticles[n].radius + zeta){
            //     roundParticles[n].forceFromCollisions[0] = (roundParticles[n].pos[0] + roundParticles[n].radius)*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
            // }
            // tmp1 = Mathf.Abs(DIM_X-1-roundParticles[n].pos[0] + roundParticles[n].radius); 
            // if(tmp1 < 2.0f*roundParticles[n].radius + zeta){
            //     roundParticles[n].forceFromCollisions[0] = -(DIM_X-1-roundParticles[n].pos[0] + roundParticles[n].radius)*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
            // }

            if(align){
                Vector2 averageMuki = new Vector2(0f,0f);
                int closeParticleCount = 0;
                for (int k = 0; k < particleCount; k++)
                {
                    if(k==n) continue;
                    for (int i = 0; i < 2; i++)
                    {
                        tmp1 = roundParticles[n].ParticleDistance(roundParticles[k]);
                        if(tmp1 < particleViewRange)
                        {
                            roundParticles[n].forceFromCollisions[i] += (roundParticles[n].pos[i] - roundParticles[k].pos[i])*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
                            averageMuki += new Vector2(Mathf.Sin(roundParticles[k].theta),Mathf.Cos(roundParticles[k].theta));
                            closeParticleCount++;
                        }
                    }
                }
                if(closeParticleCount!=0)
                {
                    averageMuki /= closeParticleCount;
                    roundParticles[n].theta = Mathf.Atan2(averageMuki[0],averageMuki[1]) + (Random.Range(-noiseAngleDeg,noiseAngleDeg)*Mathf.PI)/180f;
                }
            }
            else
            {
                for (int k = 0; k < particleCount; k++)
                {
                    if(k==n) continue;
                    for (int i = 0; i < 2; i++)
                    {
                        tmp1 = roundParticles[n].ParticleDistance(roundParticles[k]);
                        if(tmp1 < particleViewRange)
                        {
                            roundParticles[n].forceFromCollisions[i] += (roundParticles[n].pos[i] - roundParticles[k].pos[i])*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
                        }
                    }
                }
            }
            

            float wallDistance;
            float angle;
            wallDistance = Mathf.Abs(roundParticles[n].pos[1]);
            if(wallDistance < particleViewRange)
            {
                if(wallDistance < roundParticles[n].radius) roundParticles[n].pos[1] = roundParticles[n].radius;
                angle = Vector2.SignedAngle(
                    new Vector2(0,-1),
                    new Vector2(-Mathf.Sin(roundParticles[n].theta),Mathf.Cos(roundParticles[n].theta))
                );
                if(Mathf.Abs(angle) < particleFov/2f)
                {
                    float turnAngle = particleFov/2f - Mathf.Abs(angle);
                    if(angle < 0) turnAngle *= -1;
                    roundParticles[n].theta += (turnAngle*Mathf.PI)/180f;
                }
            }
            wallDistance = Mathf.Abs(DIM_Y-1-roundParticles[n].pos[1]);
            if(wallDistance < particleViewRange)
            {
                if(wallDistance < roundParticles[n].radius) roundParticles[n].pos[1] = DIM_Y-1-roundParticles[n].radius;
                angle = Vector2.SignedAngle(
                    new Vector2(0,1),
                    new Vector2(-Mathf.Sin(roundParticles[n].theta),Mathf.Cos(roundParticles[n].theta))
                );
                if(Mathf.Abs(angle) < particleFov/2f)
                {
                    float turnAngle = particleFov/2f - Mathf.Abs(angle);
                    if(angle < 0) turnAngle *= -1;
                    roundParticles[n].theta +=(turnAngle*Mathf.PI)/180f;
                }
            }
            wallDistance = Mathf.Abs(roundParticles[n].pos[0]);
            if(wallDistance < particleViewRange)
            {
                if(wallDistance < roundParticles[n].radius) roundParticles[n].pos[0] =roundParticles[n].radius;
                angle = Vector2.SignedAngle(
                    new Vector2(-1,0),
                    new Vector2(-Mathf.Sin(roundParticles[n].theta),Mathf.Cos(roundParticles[n].theta))
                );
                if(Mathf.Abs(angle) < particleFov/2f)
                {
                    float turnAngle = particleFov/2f - Mathf.Abs(angle);
                    if(angle < 0) turnAngle *= -1;
                    roundParticles[n].theta += (turnAngle*Mathf.PI)/180f;
                }
            }
            wallDistance = Mathf.Abs(DIM_X-1-roundParticles[n].pos[0]);
            if(wallDistance < particleViewRange)
            {
                if(wallDistance < roundParticles[n].radius) roundParticles[n].pos[0] = DIM_X-1-roundParticles[n].radius;
                angle = Vector2.SignedAngle(
                    new Vector2(1,0),
                    new Vector2(-Mathf.Sin(roundParticles[n].theta),Mathf.Cos(roundParticles[n].theta))
                );
                if(Mathf.Abs(angle) < particleFov/2f)
                {
                    float turnAngle = particleFov/2f - Mathf.Abs(angle);
                    if(angle < 0) turnAngle *= -1;
                    roundParticles[n].theta +=(turnAngle*Mathf.PI)/180f;
                }
            }


            

            roundParticles[n].forceFromFluid[0] = 0f;
            roundParticles[n].forceFromFluid[1] = 0f;
            roundParticles[n].torque = 0f;
            for(int m = 0; m < roundParticles[n].perimeterPointCount ; m++) 
            {
                roundParticles[n].perimeterFluidVel[m,0] = 0f;
                roundParticles[n].perimeterFluidVel[m,1] = 0f;
                angle = Vector2.Angle(
                    new Vector2(roundParticles[n].perimeterPos[m,0] - roundParticles[n].pos[0],roundParticles[n].perimeterPos[m,1] - roundParticles[n].pos[1]),
                    new Vector2(-Mathf.Sin(roundParticles[n].theta),Mathf.Cos(roundParticles[n].theta))
                );
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

                // if(m == roundParticles[n].bottomPointIndex)
                // {
                //     Vector2 vect2Center = new Vector2(
                //         roundParticles[n].pos[0] - roundParticles[n].perimeterPos[0,0],
                //         roundParticles[n].pos[1] - roundParticles[n].perimeterPos[0,1]
                //     );
                //     vect2Center.Normalize();
                //     roundParticles[n].forceOnPerimeter[m,0] -= particlePropulsion * vect2Center[0];
                //     roundParticles[n].forceOnPerimeter[m,1] -= particlePropulsion * vect2Center[1];
                // }

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

            roundParticles[n].UpdatePosVel(new float[2]{0f,gravity});
            roundParticles[n].UpdateOmegaTheta();
            roundParticles[n].UpdatePerimeter();
        }
    }
}
