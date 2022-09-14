using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LBMPoiWithPolygon : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public bool vectorFieldOn;
    public VectorField vectorField;
    public int DIM_X = 128;
    public int DIM_Y = 128;
    public float gx = 0.00001f;
    public float gy = 0f;
    public bool normalizeHeatMap = true;
    public HeatMapMode mode = HeatMapMode.Speed;

    public float tau = 0.56f;
    public float epsw = 100.0f, zeta = 1.0f;
    float tmp,tmp1,tmp2,tmp3,u2,h,nu;
    
    public float maxRho,minRho;
    public float maxSpeed,minSpeed;

    float[] cx = new float[9]{0,        1,    0,   -1,    0,     1,    -1,    -1,     1};
    float[] cy = new float[9]{0,        0,    1,    0,   -1,     1,     1,    -1,    -1};
    float[] w = new  float[9]{4f/9f,1f/9f,1f/9f,1f/9f,1f/9f,1f/36f,1f/36f,1f/36f,1f/36f};
    float[,] forceFromParticleX,forceFromParticleY;
    float[,] u,v,rho,speed;
    float[,,] f,f0,ftmp;

    public bool enableObs = true;
    public int H;
    public int W;
    public int[] obsPos = new int[2];

    public int loopCount = 1;

    public float particleRadius = 12.5f;
    public float particleDensity = 1.25f;
    public float[] particleInitPos = new float[2]{64f,64f};
    RoundParticle roundParticle;
    
    RectObstacle rectObs;

    public int lineRes = 50;
    public int modes = 4;
    public float[] modeCoeffs = new float[4];
    public float[] modeSinCoeffs = new float[4];
    public float area = 1000f;
    public float realArea;
    public GameObject dotPrefab;
    List<GameObject> dots = new List<GameObject>();
    // float dotDist = 0.5f;
    public int dotCount = 10;
    public Transform dotParent;
    public LineRenderer line;
    public Slider[] CosSliders;
    public Slider[] SinSliders;
    public Transform canvas;
    PolygonParticle polygonParticle;
    // Start is called before the first frame update
    void Start()
    {
        plotTexture = new Texture2D(DIM_X,DIM_Y);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),Vector2.zero);

        rectObs = new RectObstacle(obsPos,H,W,DIM_X);

        u = new float[DIM_X,DIM_Y];
        v = new float[DIM_X,DIM_Y];
        rho = new float[DIM_X,DIM_Y];
        speed = new float[DIM_X,DIM_Y];
        f = new float[9,DIM_X,DIM_Y];
        f0 = new float[9,DIM_X,DIM_Y];
        ftmp = new float[9,DIM_X,DIM_Y];
        forceFromParticleX = new float[DIM_X,DIM_Y];
        forceFromParticleY = new float[DIM_X,DIM_Y];

        nu = (tau - 0.5f)/3f;
        h = (float)DIM_Y-1;
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
                forceFromParticleX[i,j] = 0.0f; forceFromParticleY[i,j] = 0.0f;
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
                maxRho = Mathf.Min(maxRho,rho[i,j]);
            }
        }

        roundParticle = new RoundParticle(particleDensity,particleRadius,particleInitPos);
        DrawLine(true);
    }

    void LBMStep()
    {
        if(rectObs.h != H) rectObs.h = H;
        if(rectObs.w != W) rectObs.w = W;
        Collision();
        Streaming();
        BounceBackBoundaries();
        UpdateSpeedAndDensity();
        // ImmersedBoundary();
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
    public void OnSliderChange()
    {
        for (int i = 0; i < 3; i++)
        {
            modeCoeffs[i+1] = CosSliders[i].value;
        }
        for (int i = 0; i < 3; i++)
        {
            modeSinCoeffs[i+1] = SinSliders[i].value;
        }
        DrawLine();
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
        // roundParticle.PlotParticlePerimeter(ref plotPixels,DIM_X);
        polygonParticle.PlotParticlePerimeter(ref plotPixels,DIM_X);
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
                u2 = u[i,j]*u[i,j] + v[i,j]*v[i,j];     
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i,j] + cy[k]*v[i,j];      
                    f0[k,i,j] = w[k]*rho[i,j]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
                    f[k,i,j] += -(f[k,i,j] - f0[k,i,j])/tau;
                    // Force
                    f[k,i,j] += 3f*w[k]*rho[i,j]*(cx[k]*gx + cy[k]*gy) + 3f*w[k]*(forceFromParticleX[i,j]*cx[k] + forceFromParticleY[i,j]*cy[k]);
                } 
                // reset force from particle;
                forceFromParticleX[i,j] = 0.0f;
                forceFromParticleY[i,j] = 0.0f;
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
                    //periodic boundary condition
                    int im = (i + (int)cx[k] + DIM_X)%DIM_X; 
                    int jm = j + (int)cy[k];
                    if(jm != DIM_Y && jm !=-1) f[k,im,jm] = ftmp[k,i,j];
                } 
            }
        }
    }

    void BounceBackBoundaries()
    {
        // ftmp = (float[,,])(f.Clone());
        if(enableObs)
        {
            for (int j = rectObs.pos[1] - H/2; j <= rectObs.pos[1] + H/2; j++)
            {
                if(j<1||j>DIM_Y-1) continue;
                f[3,(rectObs.pos[0]-W/2+DIM_X)%DIM_X,j] = f[1,(rectObs.pos[0]-W/2+DIM_X)%DIM_X,j];
                f[6,(rectObs.pos[0]-W/2+DIM_X)%DIM_X,j] = f[8,(rectObs.pos[0]-W/2+DIM_X)%DIM_X,j];
                f[7,(rectObs.pos[0]-W/2+DIM_X)%DIM_X,j] = f[5,(rectObs.pos[0]-W/2+DIM_X)%DIM_X,j];

                f[1,(rectObs.pos[0]+W/2+DIM_X)%DIM_X,j] = f[3,(rectObs.pos[0]+W/2+DIM_X)%DIM_X,j];
                f[5,(rectObs.pos[0]+W/2+DIM_X)%DIM_X,j] = f[7,(rectObs.pos[0]+W/2+DIM_X)%DIM_X,j];
                f[8,(rectObs.pos[0]+W/2+DIM_X)%DIM_X,j] = f[6,(rectObs.pos[0]+W/2+DIM_X)%DIM_X,j];
            }
        }
        
        for (int i = 0; i < DIM_X; i++)
        {
            if(enableObs)
            {
                // top, bottom bounce back
                if(!rectObs.PointIsInRect(i,DIM_Y-1))
                {
                    f[4,i,DIM_Y-1] = f[2,i,DIM_Y-1];
                    f[7,i,DIM_Y-1] = f[5,i,DIM_Y-1];
                    f[8,i,DIM_Y-1] = f[6,i,DIM_Y-1]; 
                }
                if(!rectObs.PointIsInRect(i,0))
                {
                    f[2,i,0] = f[4,i,0]; 
                    f[5,i,0] = f[7,i,0]; 
                    f[6,i,0] = f[8,i,0]; 
                }

                int leftx = (rectObs.pos[0]-W/2+DIM_X)%DIM_X;
                int rightx = (rectObs.pos[0]+W/2+DIM_X)%DIM_X;
                int top = rectObs.pos[1] + H/2;
                int bottom = rectObs.pos[1] - H/2;
                bool isBoundary = false;
                if(leftx <= rightx)
                {
                    if(i>=leftx&&i<=rightx) isBoundary = true;
                }
                else
                {
                    if(i>=leftx||i<=rightx) isBoundary = true;
                }
                if(isBoundary)
                {
                    if(top>=1&&top<=DIM_Y-2)
                    {
                        f[2,i,top] = f[4,i,top]; 
                        f[5,i,top] = f[7,i,top]; 
                        f[6,i,top] = f[8,i,top]; 
                    }
                    if(bottom>=1&&bottom<=DIM_Y-2)
                    {
                        f[4,i,bottom] = f[2,i,bottom]; 
                        f[7,i,bottom] = f[5,i,bottom]; 
                        f[8,i,bottom] = f[6,i,bottom];
                    }
                    
                }
                
            }
            else
            {
                f[4,i,DIM_Y-1] = f[2,i,DIM_Y-1];
                f[7,i,DIM_Y-1] = f[5,i,DIM_Y-1];
                f[8,i,DIM_Y-1] = f[6,i,DIM_Y-1]; 
                f[2,i,0] = f[4,i,0]; 
                f[5,i,0] = f[7,i,0]; 
                f[6,i,0] = f[8,i,0]; 
            }
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

                if(enableObs)
                {
                    if(rectObs.PointIsInRect(i,j))
                    {
                        u[i,j] = 0f;
                        v[i,j] = 0f;
                    }
                }

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
        for(int n = 0; n < 1; n++) 
        { 
            roundParticle.forceFromCollisions[0] = 0f;
            roundParticle.forceFromCollisions[1] = 0f;
            tmp1 = Mathf.Abs(roundParticle.pos[1] + roundParticle.radius); 
            if(tmp1 < 2.0f*roundParticle.radius + zeta){
                roundParticle.forceFromCollisions[1] = (roundParticle.pos[1] + roundParticle.radius)*(2.0f*roundParticle.radius - tmp1 + zeta)*(2.0f*roundParticle.radius - tmp1 + zeta)/epsw;
            }
            tmp1 = Mathf.Abs(DIM_Y-roundParticle.pos[1] + roundParticle.radius); 
            if(tmp1 < 2.0f*roundParticle.radius + zeta){
                roundParticle.forceFromCollisions[1] = -(DIM_Y-roundParticle.pos[1] + roundParticle.radius)*(2.0f*roundParticle.radius - tmp1 + zeta)*(2.0f*roundParticle.radius - tmp1 + zeta)/epsw;
            }
            // tmp1 = Mathf.Abs(roundParticle.pos[0] + roundParticle.radius); 
            // if(tmp1 < 2.0f*roundParticle.radius + zeta){
            //     roundParticle.forceFromCollisions[0] = (roundParticle.pos[0] + roundParticle.radius)*(2.0f*roundParticle.radius - tmp1 + zeta)*(2.0f*roundParticle.radius - tmp1 + zeta)/epsw;
            // }
            // tmp1 = Mathf.Abs(DIM_X-1-roundParticle.pos[0] + roundParticle.radius); 
            // if(tmp1 < 2.0f*roundParticle.radius + zeta){
            //     roundParticle.forceFromCollisions[0] = -(DIM_X-1-roundParticle.pos[0] + roundParticle.radius)*(2.0f*roundParticle.radius - tmp1 + zeta)*(2.0f*roundParticle.radius - tmp1 + zeta)/epsw;
            // }

            roundParticle.forceFromFluid[0] = 0f;
            roundParticle.forceFromFluid[1] = 0f;
            roundParticle.torque = 0f;
            // 固体表面の流体の速度を計算
            for(int m = 0; m < roundParticle.perimeterPointCount ; m++) 
            {
                // uet[n,m] = 0.0f; vet[n,m] = 0.0f;
                roundParticle.perimeterFluidVel[m,0] = 0f;
                roundParticle.perimeterFluidVel[m,1] = 0f;
                for(int i = (int)roundParticle.perimeterPos[m,0] - 3; i < (int)roundParticle.perimeterPos[m,0] + 3; i++)
                {
                    for(int j = (int)roundParticle.perimeterPos[m,1] - 3; j < (int)roundParticle.perimeterPos[m,1] + 3; j++)
                    {
                        tmp1 = Mathf.Abs(roundParticle.perimeterPos[m,0] - (float)i);
                        tmp2 = Mathf.Abs(roundParticle.perimeterPos[m,1] - (float)j);
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
                        if((j<DIM_Y&&j>=0))
                        {
                            roundParticle.perimeterFluidVel[m,0] += u[(i + DIM_X)%DIM_X,j]*tmp3;
                            roundParticle.perimeterFluidVel[m,1] += v[(i + DIM_X)%DIM_X,j]*tmp3;
                        }
                    } 
                }
                roundParticle.forceOnPerimeter[m,0] = roundParticle.perimeterVel[m,0] - roundParticle.perimeterFluidVel[m,0];
                roundParticle.forceOnPerimeter[m,1] = roundParticle.perimeterVel[m,1] - roundParticle.perimeterFluidVel[m,1];

                // 固体が外部に与える力を計算
                for(int i = (int)roundParticle.perimeterPos[m,0] - 3; i < (int)roundParticle.perimeterPos[m,0] + 3; i++)
                {
                    for(int j = (int)roundParticle.perimeterPos[m,1] - 3; j < (int)roundParticle.perimeterPos[m,1] + 3; j++)
                    {
                        tmp1 = Mathf.Abs(roundParticle.perimeterPos[m,0] - (float)i);
                        tmp2 = Mathf.Abs(roundParticle.perimeterPos[m,1] - (float)j);
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
                        if((j<DIM_Y&&j>=0))
                        {
                            forceFromParticleX[(i + DIM_X)%DIM_X,j] += roundParticle.forceOnPerimeter[m,0] * tmp3 * 2.0f*Mathf.PI*roundParticle.radius/(float)roundParticle.perimeterPointCount;
                            forceFromParticleY[(i + DIM_X)%DIM_X,j] += roundParticle.forceOnPerimeter[m,1] * tmp3 * 2.0f*Mathf.PI*roundParticle.radius/(float)roundParticle.perimeterPointCount;
                        }
                    } 
                }
                roundParticle.forceFromFluid[0] += roundParticle.forceOnPerimeter[m,0];
                roundParticle.forceFromFluid[1] += roundParticle.forceOnPerimeter[m,1];
                roundParticle.torque += roundParticle.forceOnPerimeter[m,1] * (roundParticle.perimeterPos[m,0] - roundParticle.pos[0]) 
                                        - roundParticle.forceOnPerimeter[m,0] * (roundParticle.perimeterPos[m,1] - roundParticle.pos[1]);
            } 

            roundParticle.forceFromFluid[0] *= -2f*Mathf.PI*roundParticle.radius/(float)roundParticle.perimeterPointCount;  
            roundParticle.forceFromFluid[1] *= -2f*Mathf.PI*roundParticle.radius/(float)roundParticle.perimeterPointCount;  
            roundParticle.torque *= -2f*Mathf.PI*roundParticle.radius/(float)roundParticle.perimeterPointCount;  

            // roundParticle.UpdatePosVelPeriodicX(DIM_X);
            // for (int i = 0; i < 2; i++)
            // {
            //     roundParticle.vel[i] = (1f + 1f/roundParticle.density) * roundParticle.prevVel1[i]
            //                             - 1f/roundParticle.density * roundParticle.prevVel2[i]
            //                             + (roundParticle.forceFromFluid[i] + roundParticle.forceFromCollisions[i])/roundParticle.mass;
            //                             // + (1f - 1f/roundParticle.density) * gravity[i];
            //     roundParticle.pos[i] += (roundParticle.vel[i] + roundParticle.prevVel1[i])/2f;
            //     roundParticle.prevVel2[i] = roundParticle.prevVel1[i];
            //     roundParticle.prevVel1[i] = roundParticle.vel[i];
            // }
            // roundParticle.omega = (1f + 1f/roundParticle.density) * roundParticle.prevOmega1 
            //                     - 1f/roundParticle.density * roundParticle.prevOmega2
            //                     + roundParticle.torque/roundParticle.momentOfInertia;
            // roundParticle.theta += (roundParticle.omega - roundParticle.prevOmega1)/2f;

            // roundParticle.prevOmega2 = roundParticle.prevOmega1;
            // roundParticle.prevOmega1 = roundParticle.omega;
            // roundParticle.UpdatePerimeterPeriodicX(DIM_X);
        }
    }

    void DrawLine(bool initParticle = false)
    {
        List<Vector3> points = new List<Vector3>();
        float circumferenceProgressPerStep = (float)1/lineRes;
        float TAU = 2*Mathf.PI;
        float radianProgressPerStep = circumferenceProgressPerStep*TAU;
        float scale;
        float sqrdSum = 0f;
        for (int i = 0; i < lineRes; i++)
        {
            float currentRadian = radianProgressPerStep*i;
            float newRadius = 0f;
            for (int j = 0; j < modes; j++)
            {
                newRadius += modeCoeffs[j]*Mathf.Cos(j*(currentRadian)) + modeSinCoeffs[j]*Mathf.Sin(j*(currentRadian));
            }
            sqrdSum += newRadius * newRadius;
        }
        scale = Mathf.Sqrt((2*area)/(sqrdSum * lineRes * Mathf.Sin(TAU/lineRes)));

        realArea = 0f;
        for(int i = 0; i < lineRes; i++)
        {
            float currentRadian = radianProgressPerStep*i;
            Vector3 newPoint = new Vector3(Mathf.Cos(currentRadian), Mathf.Sin(currentRadian),0);
            float newRadius = 0f;
            for (int j = 0; j < modes; j++)
            {
                newRadius += modeCoeffs[j]*Mathf.Cos(j*(currentRadian)) + modeSinCoeffs[j]*Mathf.Sin(j*(currentRadian));
            }
            newRadius *= scale;
            newPoint *= newRadius;
            realArea += newRadius*newRadius;
            points.Add(newPoint);
        }
        realArea *= (lineRes * Mathf.Sin(TAU/lineRes))/2f;

        float perimeterLength = 2f * Mathf.PI * modeCoeffs[0] * scale;
        float segmentLength = perimeterLength/dotCount;
        float theta = 0f;
        float dtheta,r;
        int spawnedDots = 0;
        Vector3 spawnPos;
        while (theta <= 2*Mathf.PI)
        {
            r = 0f;
            for (int i = 0; i < modes; i++)
            {
                r += modeCoeffs[i]*Mathf.Cos(i*theta) + modeSinCoeffs[i]*Mathf.Sin(i*theta);
            }
            r *= scale;
            if(r <= Mathf.Epsilon)
            {
                theta += 0.1f * Mathf.PI;
                continue;
            }
            spawnPos = new Vector3(Mathf.Cos(theta), Mathf.Sin(theta),0);
            spawnPos *= r;
            if(dots.Count < spawnedDots + 1)
            {
                dots.Add(
                    Instantiate(dotPrefab,spawnPos,Quaternion.identity,dotParent)
                );
            }
            else
            {
                dots[spawnedDots].transform.position = spawnPos;
            }
            spawnedDots++;

            dtheta = segmentLength/r;
            theta += dtheta;
        }
        if(spawnedDots<dots.Count)
        {
            for (int i = 0; i < dots.Count - spawnedDots; i++)
            {
                Destroy(dots[dots.Count-1]);
                dots.RemoveAt(dots.Count-1);
            }
        }
        int pointCount = dots.Count;
        float[,] perimeterPos = new float[pointCount,2];
        for (int i = 0; i < pointCount; i++)
        {
            perimeterPos[i,0] = dots[i].transform.position.x/canvas.localScale.x + plotImage.transform.GetComponent<RectTransform>().sizeDelta.x*plotImage.transform.localScale.x/2f;
            perimeterPos[i,1] = dots[i].transform.position.y/canvas.localScale.y + plotImage.transform.GetComponent<RectTransform>().sizeDelta.y*plotImage.transform.localScale.y/2f;
            perimeterPos[i,0] *= DIM_X/(plotImage.transform.GetComponent<RectTransform>().sizeDelta.x*plotImage.transform.localScale.x);
            perimeterPos[i,1] *= DIM_Y/(plotImage.transform.GetComponent<RectTransform>().sizeDelta.y*plotImage.transform.localScale.y);
        }
        if(initParticle)polygonParticle = new PolygonParticle(pointCount,perimeterPos);
        else {
            polygonParticle.perimeterPointCount = pointCount;
            polygonParticle.perimeterPos = perimeterPos;
        }
        points.Add(points[0]);
        line.positionCount = points.Count;
        line.SetPositions(points.ToArray());
    }
}

