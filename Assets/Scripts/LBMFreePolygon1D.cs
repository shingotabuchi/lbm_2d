using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LBMFreePolygon1D : MonoBehaviour
{
    public Image plotImage;
    public Image tracePlotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    Texture2D tracePlotTexture;
    Color[] tracePlotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public bool vectorFieldOn;
    public VectorField vectorField;
    public int DIM_X = 128;
    public int DIM_Y = 128;
    public float gx = 0f;
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
    RectObstacle rectObs;

    public int lineRes = 50;
    public int modes = 4;
    public float[] modeCoeffs = new float[4];
    public float[] modeSinCoeffs = new float[4];
    public float area = 1000f;
    public float realArea;
    public int dotCount = 10;
    public LineRenderer line;
    public Slider[] CosSliders;
    public Slider[] SinSliders;
    public Transform canvas;
    public Color perimColor;
    public Color traceColor;

    PolygonParticle polygonParticle;

    List<Vector3> dotPositions = new List<Vector3>();

    public float pTheta;

    public bool falling = false;
    // Start is called before the first frame update
    void Start()
    {
        plotTexture = new Texture2D(DIM_X,DIM_Y);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),Vector2.zero);

        tracePlotTexture = new Texture2D(DIM_X,DIM_Y);
        tracePlotTexture.filterMode = FilterMode.Point;
        tracePlotPixels = tracePlotTexture.GetPixels();
        for (int i = 0; i < tracePlotPixels.Length; i++)
        {
            tracePlotPixels[i] = new Color(0,0,0,0);
        }
        tracePlotTexture.SetPixels(tracePlotPixels);
        tracePlotTexture.Apply();
        tracePlotImage.sprite = Sprite.Create(tracePlotTexture, new Rect(0,0,DIM_X,DIM_Y),Vector2.zero);

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

        // polygonParticle = new RoundParticle(particleDensity,particleRadius,particleInitPos);
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
        // pTheta = polygonParticle.theta;
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
        // DrawLine();
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
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                tracePlotPixels[(int)polygonParticle.pos[0]+i + (int)(polygonParticle.pos[1]+j)*tracePlotTexture.width] = traceColor;
            }
        }
        tracePlotTexture.SetPixels(tracePlotPixels);
        tracePlotTexture.Apply();
        // polygonParticle.PlotParticlePerimeter(ref plotPixels,DIM_X);
        polygonParticle.PlotParticlePerimeter(ref plotPixels,DIM_X,perimColor);
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
                    if(falling) f[k,i,j] += 3f*w[k]*(forceFromParticleX[i,j]*cx[k] + forceFromParticleY[i,j]*cy[k]);
                    else f[k,i,j] += 3f*w[k]*rho[i,j]*(cx[k]*gx + cy[k]*gy) + 3f*w[k]*(forceFromParticleX[i,j]*cx[k] + forceFromParticleY[i,j]*cy[k]);
                } 
                // reset force from particle;
                forceFromParticleX[i,j] = 0.0f;
                forceFromParticleY[i,j] = 0.0f;
            }   
        }
    }

    void Streaming()
    {
        int im,jm;
        ftmp = (float[,,])(f.Clone());
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            { 
                for(int k = 0; k <= 8; k++)
                {
                    if(falling)
                    {
                        im = i + (int)cx[k]; 
                        jm = j + (int)cy[k];
                        if((jm!=DIM_Y&&jm!=-1) && (im!=DIM_X&&im!=-1))
                        {
                            f[k,im,jm] = ftmp[k,i,j];
                        }
                    }
                    else
                    {
                        //periodic boundary condition
                        im = (i + (int)cx[k] + DIM_X)%DIM_X; 
                        jm = j + (int)cy[k];
                        if(jm != DIM_Y && jm !=-1) f[k,im,jm] = ftmp[k,i,j];
                    }
                } 
            }
        }
    }

    void BounceBackBoundaries()
    {
        // ftmp = (float[,,])(f.Clone());
        if(falling)
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
            return;
        }
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
        polygonParticle.forceFromFluid[0] = 0f;
        polygonParticle.forceFromFluid[1] = 0f;
        polygonParticle.torque = 0f;
        // 固体表面の流体の速度を計算
        for(int m = 0; m < polygonParticle.perimeterPointCount ; m++) 
        {
            // uet[n,m] = 0.0f; vet[n,m] = 0.0f;
            polygonParticle.perimeterFluidVel[m,0] = 0f;
            polygonParticle.perimeterFluidVel[m,1] = 0f;
            for(int i = (int)polygonParticle.perimeterPos[m,0] - 3; i < (int)polygonParticle.perimeterPos[m,0] + 3; i++)
            {
                for(int j = (int)polygonParticle.perimeterPos[m,1] - 3; j < (int)polygonParticle.perimeterPos[m,1] + 3; j++)
                {
                    tmp1 = Mathf.Abs(polygonParticle.perimeterPos[m,0] - (float)i);
                    tmp2 = Mathf.Abs(polygonParticle.perimeterPos[m,1] - (float)j);
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
                        polygonParticle.perimeterFluidVel[m,0] += u[(i + DIM_X)%DIM_X,j]*tmp3;
                        polygonParticle.perimeterFluidVel[m,1] += v[(i + DIM_X)%DIM_X,j]*tmp3;
                    }
                } 
            }
            polygonParticle.forceOnPerimeter[m,0] = polygonParticle.perimeterVel[m,0] - polygonParticle.perimeterFluidVel[m,0];
            polygonParticle.forceOnPerimeter[m,1] = polygonParticle.perimeterVel[m,1] - polygonParticle.perimeterFluidVel[m,1];

            // 固体が外部に与える力を計算
            for(int i = (int)polygonParticle.perimeterPos[m,0] - 3; i < (int)polygonParticle.perimeterPos[m,0] + 3; i++)
            {
                for(int j = (int)polygonParticle.perimeterPos[m,1] - 3; j < (int)polygonParticle.perimeterPos[m,1] + 3; j++)
                {
                    tmp1 = Mathf.Abs(polygonParticle.perimeterPos[m,0] - (float)i);
                    tmp2 = Mathf.Abs(polygonParticle.perimeterPos[m,1] - (float)j);
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
                        forceFromParticleX[(i + DIM_X)%DIM_X,j] += polygonParticle.forceOnPerimeter[m,0] * tmp3 * 0.5f;
                        forceFromParticleY[(i + DIM_X)%DIM_X,j] += polygonParticle.forceOnPerimeter[m,1] * tmp3 * 0.5f;
                    }
                } 
            }
            polygonParticle.forceFromFluid[0] += polygonParticle.forceOnPerimeter[m,0];
            polygonParticle.forceFromFluid[1] += polygonParticle.forceOnPerimeter[m,1];
            polygonParticle.torque += polygonParticle.forceOnPerimeter[m,1] * (polygonParticle.perimeterPos[m,0] - polygonParticle.pos[0]) 
                                    - polygonParticle.forceOnPerimeter[m,0] * (polygonParticle.perimeterPos[m,1] - polygonParticle.pos[1]);
        } 

        polygonParticle.forceFromFluid[0] *= -0.5f;  
        polygonParticle.forceFromFluid[1] *= -0.5f; 
        polygonParticle.torque *= -0.5f; 

        if(falling) polygonParticle.UpdatePosVel(new float[2]{gx,gy});
        else polygonParticle.UpdatePosVel();
        polygonParticle.UpdateOmegaTheta();
        polygonParticle.theta = pTheta * Mathf.PI;
        UpdatePolygonPerimeter();
    }
    void UpdatePolygonPerimeter()
    {
        float circumferenceProgressPerStep = (float)1/lineRes;
        float TAU = 2*Mathf.PI;
        float radianProgressPerStep = circumferenceProgressPerStep*TAU;
        float scale;
        float sqrdSum = 0f;
        float polygonTheta = 0f;
        polygonTheta = -polygonParticle.theta;
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

        float perimeterLength = 2f * Mathf.PI * modeCoeffs[0] * scale;
        float segmentLength = 1f/(2f* (DIM_X/((canvas.localScale.x)*(plotImage.transform.GetComponent<RectTransform>().sizeDelta.x*plotImage.transform.localScale.x))));
        dotCount = (int)(perimeterLength/segmentLength);
        float theta = 0f;
        float dtheta,r;
        int spawnedDots = 0;
        Vector3 spawnPos;
        while (theta <= 2*Mathf.PI)
        {
            r = 0f;
            for (int i = 0; i < modes; i++)
            {
                r += modeCoeffs[i]*Mathf.Cos(i*(theta)) + modeSinCoeffs[i]*Mathf.Sin(i*(theta));
            }
            r *= scale;
            if(r <= Mathf.Epsilon)
            {
                theta += 0.1f * Mathf.PI;
                continue;
            }
            spawnPos = new Vector3(Mathf.Cos((theta+polygonTheta)), Mathf.Sin((theta+polygonTheta)),0);
            spawnPos *= r;
            if(dotPositions.Count < spawnedDots + 1)
            {
                dotPositions.Add(spawnPos);
            }
            else
            {
                dotPositions[spawnedDots] = spawnPos;
            }
            spawnedDots++;

            dtheta = segmentLength/r;
            theta += dtheta;
        }
        if(spawnedDots<dotPositions.Count)
        {
            for (int i = 0; i < dotPositions.Count - spawnedDots; i++)
            {
                dotPositions.RemoveAt(dotPositions.Count-1);
            }
        }
        polygonParticle.perimeterPointCount =  dotPositions.Count;
        polygonParticle.perimeterPos = new float[polygonParticle.perimeterPointCount,2];
        for (int i = 0; i < polygonParticle.perimeterPointCount; i++)
        {
            float newposx = dotPositions[i].x/canvas.localScale.x;
            newposx *= DIM_X/(plotImage.transform.GetComponent<RectTransform>().sizeDelta.x*plotImage.transform.localScale.x);
            float newposy = dotPositions[i].y/canvas.localScale.y;
            newposy *= DIM_Y/(plotImage.transform.GetComponent<RectTransform>().sizeDelta.y*plotImage.transform.localScale.y);
            polygonParticle.perimeterPos[i,0] = polygonParticle.pos[0] + newposx;
            polygonParticle.perimeterPos[i,1] = polygonParticle.pos[1] + newposy;
            polygonParticle.perimeterVel[i,0] = polygonParticle.vel[0] - polygonParticle.omega*(polygonParticle.perimeterPos[i,1] - polygonParticle.pos[1]);
            polygonParticle.perimeterVel[i,1] = polygonParticle.vel[1] + polygonParticle.omega*(polygonParticle.perimeterPos[i,0] - polygonParticle.pos[0]);
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
        float polygonTheta = 0f;
        if(initParticle) polygonTheta = 0f;
        else polygonTheta = -polygonParticle.theta;
        for (int i = 0; i < lineRes; i++)
        {
            float currentRadian = radianProgressPerStep*i + polygonTheta;
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
            float currentRadian = radianProgressPerStep*i + polygonTheta;
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
        // float segmentLength = perimeterLength/dotCount;
        float segmentLength = 1f/(2f* (DIM_X/((canvas.localScale.x)*(plotImage.transform.GetComponent<RectTransform>().sizeDelta.x*plotImage.transform.localScale.x))));
        dotCount = (int)(perimeterLength/segmentLength);
        float theta = 0f;
        float dtheta,r;
        int spawnedDots = 0;
        Vector3 spawnPos;
        while (theta <= 2*Mathf.PI)
        {
            r = 0f;
            for (int i = 0; i < modes; i++)
            {
                r += modeCoeffs[i]*Mathf.Cos(i*(theta+polygonTheta)) + modeSinCoeffs[i]*Mathf.Sin(i*(theta+polygonTheta));
            }
            r *= scale;
            if(r <= Mathf.Epsilon)
            {
                theta += 0.1f * Mathf.PI;
                continue;
            }
            spawnPos = new Vector3(Mathf.Cos((theta+polygonTheta)), Mathf.Sin((theta+polygonTheta)),0);
            spawnPos *= r;


            // if(dots.Count < spawnedDots + 1)
            if(dotPositions.Count < spawnedDots + 1)
            {
                dotPositions.Add(spawnPos);
                // dots.Add(
                //     Instantiate(dotPrefab,spawnPos,Quaternion.identity,dotParent)
                // );
            }
            else
            {
                // dots[spawnedDots].transform.position = spawnPos;
                dotPositions[spawnedDots] = spawnPos;
            }
            spawnedDots++;

            dtheta = segmentLength/r;
            theta += dtheta;
        }
        // if(spawnedDots<dots.Count)
        if(spawnedDots<dotPositions.Count)
        {
            for (int i = 0; i < dotPositions.Count - spawnedDots; i++)
            {
                // Destroy(dots[dots.Count-1]);
                // dots.RemoveAt(dots.Count-1);
                dotPositions.RemoveAt(dotPositions.Count-1);
            }
        }
        // int pointCount = dots.Count;
        if(initParticle){
            // polygonParticle = new PolygonParticle(spawnedDots,particleDensity,area,new float[2]{(float)(DIM_X)/2f,(float)(DIM_Y)/2f});
            polygonParticle = new PolygonParticle(spawnedDots,particleDensity,area,new float[2]{(float)(DIM_X)/2f,350});
        }
        polygonParticle.perimeterPointCount =  dotPositions.Count;
        polygonParticle.perimeterPos = new float[polygonParticle.perimeterPointCount,2];
        for (int i = 0; i < polygonParticle.perimeterPointCount; i++)
        {
            // perimeterPos[i,0] = dots[i].transform.position.x/canvas.localScale.x + plotImage.transform.GetComponent<RectTransform>().sizeDelta.x*plotImage.transform.localScale.x/2f;
            // perimeterPos[i,1] = dots[i].transform.position.y/canvas.localScale.y + plotImage.transform.GetComponent<RectTransform>().sizeDelta.y*plotImage.transform.localScale.y/2f;
            float newposx = dotPositions[i].x/canvas.localScale.x;
            newposx *= DIM_X/(plotImage.transform.GetComponent<RectTransform>().sizeDelta.x*plotImage.transform.localScale.x);
            float newposy = dotPositions[i].y/canvas.localScale.y;
            newposy *= DIM_Y/(plotImage.transform.GetComponent<RectTransform>().sizeDelta.y*plotImage.transform.localScale.y);
            polygonParticle.perimeterPos[i,0] = polygonParticle.pos[0] + newposx;
            polygonParticle.perimeterPos[i,1] = polygonParticle.pos[1] + newposy;
        }
        // print(polygonParticle.pos[0]);
        // print(polygonParticle.pos[1]);
        points.Add(points[0]);
        line.positionCount = points.Count;
        line.SetPositions(points.ToArray());
        line.transform.gameObject.SetActive(false);
    }
}

