using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Lbmpoi : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public VectorField vectorField;
    public int DIM = 64;
    public float gx = 0.00001f;
    public float gy = 0f;
    public HeatMapMode mode = HeatMapMode.Speed;

    public float tau = 0.56f;
    float tmp,u2,h,nu,norm,nb,maxrho,prevMaxSpeed;
    public float maxSpeed;
    float[] cx = new float[9]{0,        1,    0,   -1,    0,     1,    -1,    -1,     1};
    float[] cy = new float[9]{0,        0,    1,    0,   -1,     1,     1,    -1,    -1};
    float[] w = new  float[9]{4f/9f,1f/9f,1f/9f,1f/9f,1f/9f,1f/36f,1f/36f,1f/36f,1f/36f};
    float[,] u,v,rho,speed;
    float[,,] f,f0,ftmp;

    public bool enableObs = true;
    public int H;
    public int W;
    public int[] obsPos = new int[2];

    public int loopCount = 1;
    
    RectObstacle rectObs;
    // Start is called before the first frame update
    void Start()
    {
        plotTexture = new Texture2D(DIM,DIM);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM,DIM),Vector2.zero);

        rectObs = new RectObstacle(obsPos,H,W,DIM);

        u = new float[DIM,DIM];
        v = new float[DIM,DIM];
        rho = new float[DIM,DIM];
        speed = new float[DIM,DIM];
        f = new float[9,DIM,DIM];
        f0 = new float[9,DIM,DIM];
        ftmp = new float[9,DIM,DIM];

        nu = (tau - 0.5f)/3f;
        h = (float)DIM-1;
        maxSpeed = 0f;
        maxrho = 0f;
        for(int i = 0; i < DIM; i++)
        { 
            for(int j = 0; j < DIM; j++)
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
                maxrho = Mathf.Max(maxrho,rho[i,j]);
            }
        }
        // float diffsqrd = 1f;
        // int t = 0;
        // while(diffsqrd > 1e-14 && t<10000)
        // {
        //     t++;
        //     prevMaxSpeed = maxSpeed;
        //     LBMStep();
        //     diffsqrd = (maxSpeed - prevMaxSpeed)*(maxSpeed - prevMaxSpeed);
        // }
        // print(t);
    }

    void LBMStep()
    {
        if(rectObs.h != H) rectObs.h = H;
        if(rectObs.w != W) rectObs.w = W;
        Collision();
        Streaming();
        BounceBackBoundaries();
        UpdateSpeedAndDensity();
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
            if(mode == HeatMapMode.Speed)
            {
                plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM,i/DIM],maxSpeed);
            }
            else if(mode == HeatMapMode.Density)
            {
                plotPixels[i] = colorHeatMap.GetColorForValue(rho[i%DIM,i/DIM],maxrho);
            }
        }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
        vectorField.UpdateVectors(u,v,DIM,maxSpeed);
    }

    void Collision()
    {
        for(int i = 0; i < DIM; i++)
        { 
            for(int j = 0; j < DIM; j++)
            {
                u2 = u[i,j]*u[i,j] + v[i,j]*v[i,j];     
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i,j] + cy[k]*v[i,j];      
                    f0[k,i,j] = w[k]*rho[i,j]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
                    f[k,i,j] += -(f[k,i,j] - f0[k,i,j])/tau;
                    // Force
                    f[k,i,j] += 3f*w[k]*rho[i,j]*(cx[k]*gx + cy[k]*gy);
                } 
            }   
        }
    }

    void Streaming()
    {
        ftmp = (float[,,])(f.Clone());
        for(int i = 0; i < DIM; i++)
        { 
            for(int j = 0; j < DIM; j++)
            { 
                for(int k = 0; k <= 8; k++)
                {
                    //periodic boundary condition
                    int im = (i + (int)cx[k] + DIM)%DIM; 
                    int jm = j + (int)cy[k];
                    if(jm != DIM && jm !=-1) f[k,im,jm] = ftmp[k,i,j];
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
                if(j<1||j>DIM-1) continue;
                f[3,(rectObs.pos[0]-W/2+DIM)%DIM,j] = f[1,(rectObs.pos[0]-W/2+DIM)%DIM,j];
                f[6,(rectObs.pos[0]-W/2+DIM)%DIM,j] = f[8,(rectObs.pos[0]-W/2+DIM)%DIM,j];
                f[7,(rectObs.pos[0]-W/2+DIM)%DIM,j] = f[5,(rectObs.pos[0]-W/2+DIM)%DIM,j];

                f[1,(rectObs.pos[0]+W/2+DIM)%DIM,j] = f[3,(rectObs.pos[0]+W/2+DIM)%DIM,j];
                f[5,(rectObs.pos[0]+W/2+DIM)%DIM,j] = f[7,(rectObs.pos[0]+W/2+DIM)%DIM,j];
                f[8,(rectObs.pos[0]+W/2+DIM)%DIM,j] = f[6,(rectObs.pos[0]+W/2+DIM)%DIM,j];
            }
        }
        
        for (int i = 0; i < DIM; i++)
        {
            if(enableObs)
            {
                // top, bottom bounce back
                if(!rectObs.PointIsInRect(i,DIM-1))
                {
                    f[4,i,DIM-1] = f[2,i,DIM-1];
                    f[7,i,DIM-1] = f[5,i,DIM-1];
                    f[8,i,DIM-1] = f[6,i,DIM-1]; 
                }
                if(!rectObs.PointIsInRect(i,0))
                {
                    f[2,i,0] = f[4,i,0]; 
                    f[5,i,0] = f[7,i,0]; 
                    f[6,i,0] = f[8,i,0]; 
                }

                int leftx = (rectObs.pos[0]-W/2+DIM)%DIM;
                int rightx = (rectObs.pos[0]+W/2+DIM)%DIM;
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
                    if(top>=1&&top<=DIM-2)
                    {
                        f[2,i,top] = f[4,i,top]; 
                        f[5,i,top] = f[7,i,top]; 
                        f[6,i,top] = f[8,i,top]; 
                    }
                    if(bottom>=1&&bottom<=DIM-2)
                    {
                        f[4,i,bottom] = f[2,i,bottom]; 
                        f[7,i,bottom] = f[5,i,bottom]; 
                        f[8,i,bottom] = f[6,i,bottom];
                    }
                    
                }
                
            }
            else
            {
                f[4,i,DIM-1] = f[2,i,DIM-1];
                f[7,i,DIM-1] = f[5,i,DIM-1];
                f[8,i,DIM-1] = f[6,i,DIM-1]; 
                f[2,i,0] = f[4,i,0]; 
                f[5,i,0] = f[7,i,0]; 
                f[6,i,0] = f[8,i,0]; 
            }
        }
    }

    void UpdateSpeedAndDensity()
    {
        maxSpeed = 0f;
        maxrho = 0f;
        for(int i = 0; i < DIM; i++)
        { 
            for(int j = 0; j < DIM; j++)
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
                maxrho = Mathf.Max(maxrho,rho[i,j]);
            } 
        }
    }
}
