using System.Collections;
using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class IBLBMSingle : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public VectorField vectorField;
    public bool vectorFieldOn;
    public HeatMapMode mode = HeatMapMode.Speed;
    public bool normalizeHeatMap;
    public int DIM_X = 100;
    public int DIM_Y = 400;
    public int particleCount = 1;

    public float maxRho,minRho;
    public float maxSpeed,minSpeed;

    float[] cx = new float[9]{0, 1,    0,   -1,    0,     1,    -1,    -1,     1};
    float[] cy = new float[9]{0, 0,    1,    0,   -1,     1,     1,    -1,    -1};
    float[] w = new float[9]{4f/9f,1f/9f,1f/9f,1f/9f,1f/9f,1f/36f,1f/36f,1f/36f,1f/36f};
    float[,] rho, u, v, speed, fx ,fy;
    float[,,] f, f0, ftmp;
    float tmp,u2,nu,tmp1, tmp2, tmp3,dt;
    public float particleDensity = 1.25f;
    public float particleRadius = 12f;
    public float epsw = 100.0f, zeta = 1.0f, tau = 0.53f;

    public int loopCount = 1;
    public float gRate = 6f;
    float[] gravity = new float[2];
    RoundParticle roundParticle;

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
        f = new float[9,DIM_X,DIM_Y];
        f0 = new float[9,DIM_X,DIM_Y];
        ftmp = new float[9,DIM_X,DIM_Y];

        nu = (tau - 0.5f)/3.0f;

        dt = nu/((float)(DIM_X-1)/2.0f)/((float)(DIM_X-1)/2.0f)/0.1f;
        gravity[1] = -gRate*981.0f*(float)(DIM_X-1)*dt*dt;
        gravity[0] = 0.0f;

        roundParticle = new RoundParticle(particleDensity,particleRadius,new float[2]{50f,300f});
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
        BouncebackBoundaries();
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
        roundParticle.PlotParticlePerimeter(ref plotPixels,DIM_X);
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
                    // // periodic boundary
                    // int im = (i + (int)cx[k] + DIM_X)%DIM_X; 
                    // int jm = (j + (int)cy[k] + DIM_Y)%DIM_Y;
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
        for(int n = 0; n < 1; n++) 
        { 
            roundParticle.forceFromCollisions[0] = 0f;
            roundParticle.forceFromCollisions[1] = 0f;
            tmp1 = Mathf.Abs(roundParticle.pos[1] + roundParticle.radius); 
            if(tmp1 < 2.0f*roundParticle.radius + zeta){
                roundParticle.forceFromCollisions[1] = (roundParticle.pos[1] + roundParticle.radius)*(2.0f*roundParticle.radius - tmp1 + zeta)*(2.0f*roundParticle.radius - tmp1 + zeta)/epsw;
            }
            tmp1 = Mathf.Abs(DIM_Y-1-roundParticle.pos[1] + roundParticle.radius); 
            if(tmp1 < 2.0f*roundParticle.radius + zeta){
                roundParticle.forceFromCollisions[1] = -(DIM_Y-1-roundParticle.pos[1] + roundParticle.radius)*(2.0f*roundParticle.radius - tmp1 + zeta)*(2.0f*roundParticle.radius - tmp1 + zeta)/epsw;
            }
            tmp1 = Mathf.Abs(roundParticle.pos[0] + roundParticle.radius); 
            if(tmp1 < 2.0f*roundParticle.radius + zeta){
                roundParticle.forceFromCollisions[0] = (roundParticle.pos[0] + roundParticle.radius)*(2.0f*roundParticle.radius - tmp1 + zeta)*(2.0f*roundParticle.radius - tmp1 + zeta)/epsw;
            }
            tmp1 = Mathf.Abs(DIM_X-1-roundParticle.pos[0] + roundParticle.radius); 
            if(tmp1 < 2.0f*roundParticle.radius + zeta){
                roundParticle.forceFromCollisions[0] = -(DIM_X-1-roundParticle.pos[0] + roundParticle.radius)*(2.0f*roundParticle.radius - tmp1 + zeta)*(2.0f*roundParticle.radius - tmp1 + zeta)/epsw;
            }

            // fxc[n] = 0.0f;
            // fyc[n] = 0.0f;
            // if(tmp1 < 2.0f*rp + zeta){
            //     fyc[n] += (yp[n] + rp)*(2.0f*rp - tmp1 + zeta)*(2.0f*rp - tmp1 + zeta)/epsw;
            // }

            // fxp[n] = 0.0f; fyp[n] = 0.0f; torq[n] = 0.0f;
            roundParticle.forceFromFluid[0] = 0f;
            roundParticle.forceFromFluid[1] = 0f;
            roundParticle.torque = 0f;
            for(int m = 0; m < roundParticle.perimeterPointCount ; m++) 
            {
                // uet[n,m] = 0.0f; vet[n,m] = 0.0f;
                roundParticle.perimeterFluidVel[m,0] = 0f;
                roundParticle.perimeterFluidVel[m,1] = 0f;
                // 固体表面の速度を計算
                // for(int i = (int)xe[n,m] - 3; i < (int)xe[n,m] + 3; i++)
                // {
                //     for(int j = (int)ye[n,m] - 3; j < (int)ye[n,m] + 3; j++)
                //     {
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
                        if((j<DIM_Y&&j>=0) && (i<DIM_X&&i>=0))
                        {
                            // uet[n,m] += u[i,j]*tmp3;
                            // vet[n,m] += v[i,j]*tmp3;
                            roundParticle.perimeterFluidVel[m,0] += u[i,j]*tmp3;
                            roundParticle.perimeterFluidVel[m,1] += v[i,j]*tmp3;
                        }
                    } 
                }
                // fxe[n,m] = ue[n,m] - uet[n,m];
                // fye[n,m] = ve[n,m] - vet[n,m];
                roundParticle.forceOnPerimeter[m,0] = roundParticle.perimeterVel[m,0] - roundParticle.perimeterFluidVel[m,0];
                roundParticle.forceOnPerimeter[m,1] = roundParticle.perimeterVel[m,1] - roundParticle.perimeterFluidVel[m,1];

                // 固体が外部に与える力を計算
                // for(int i = (int)xe[n,m] - 3; i < (int)xe[n,m] + 3; i++)
                // {
                //     for(int j = (int)ye[n,m] - 3; j < (int)ye[n,m] + 3; j++)
                //     {
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
                        if((j<DIM_Y&&j>=0) && (i<DIM_X&&i>=0))
                        {
                            fx[i,j] += roundParticle.forceOnPerimeter[m,0] * tmp3 * 2.0f*Mathf.PI*roundParticle.radius/(float)roundParticle.perimeterPointCount;
                            fy[i,j] += roundParticle.forceOnPerimeter[m,1] * tmp3 * 2.0f*Mathf.PI*roundParticle.radius/(float)roundParticle.perimeterPointCount;
                        }
                    } 
                }
                roundParticle.forceFromFluid[0] += roundParticle.forceOnPerimeter[m,0];
                roundParticle.forceFromFluid[1] += roundParticle.forceOnPerimeter[m,1];
                roundParticle.torque += roundParticle.forceOnPerimeter[m,1] * (roundParticle.perimeterPos[m,0] - roundParticle.pos[0]) 
                                        - roundParticle.forceOnPerimeter[m,0] * (roundParticle.perimeterPos[m,1] - roundParticle.pos[1]);
                // torq[n] = torq[n] + fye[n,m]*(xe[n,m] - xp[n]) - fxe[n,m]*(ye[n,m] - yp[n]);
            } 

            // fxp[n] =  -fxp[n]*2.0f*Mathf.PI*rp/(float)ne;
            // fyp[n] =  -fyp[n]*2.0f*Mathf.PI*rp/(float)ne;
            // torq[n] = -torq[n]*2.0f*Mathf.PI*rp/(float)ne;

            roundParticle.forceFromFluid[0] *= -2f*Mathf.PI*roundParticle.radius/(float)roundParticle.perimeterPointCount;  
            roundParticle.forceFromFluid[1] *= -2f*Mathf.PI*roundParticle.radius/(float)roundParticle.perimeterPointCount;  
            roundParticle.torque *= -2f*Mathf.PI*roundParticle.radius/(float)roundParticle.perimeterPointCount;  

            for (int i = 0; i < 2; i++)
            {
                roundParticle.vel[i] = (1f + 1f/roundParticle.density) * roundParticle.prevVel1[i]
                                        - 1f/roundParticle.density * roundParticle.prevVel2[i]
                                        + (roundParticle.forceFromFluid[i] + roundParticle.forceFromCollisions[i])/roundParticle.mass
                                        + (1f - 1f/roundParticle.density) * gravity[i];
                roundParticle.pos[i] += (roundParticle.vel[i] + roundParticle.prevVel1[i])/2f;
                roundParticle.prevVel2[i] = roundParticle.prevVel1[i];
                roundParticle.prevVel1[i] = roundParticle.vel[i];
            }
            
            // up[n]  = (1.0f + 1.0f/rhop)*up1[n]    - 1.0f/rhop*up2[n]
            //         + (fxp[n] + fxc[n])/rhop/mp + (1.0f - 1.0f/rhop)*gx;
            // vp[n]  = (1.0f + 1.0f/rhop)*vp1[n]    - 1.0f/rhop*vp2[n]
            //         + (fyp[n] + fyc[n])/rhop/mp + (1.0f - 1.0f/rhop)*gy;
            // omega[n]= (1.0f + 1.0f/rhop)*omega1[n] - 1.0f/rhop*omega2[n]
            //         + torq[n]/rhop/iip;

            roundParticle.omega = (1f + 1f/roundParticle.density) * roundParticle.prevOmega1 
                                - 1f/roundParticle.density * roundParticle.prevOmega2
                                + roundParticle.torque/roundParticle.momentOfInertia;
            roundParticle.theta += (roundParticle.omega - roundParticle.prevOmega1)/2f;

            roundParticle.prevOmega2 = roundParticle.prevOmega1;
            roundParticle.prevOmega1 = roundParticle.omega;
                // xp[n]=    xp[n] + (   up[n] +    up1[n])*0.5f;
                // yp[n]=    yp[n] + (   vp[n] +    vp1[n])*0.5f;
            // theta[n]= theta[n] + (omega[n] + omega1[n])*0.5f;
                // up2[n] =    up1[n];    up1[n] =    up[n];
                // vp2[n] =    vp1[n];    vp1[n] =    vp[n];
            // omega2[n] = omega1[n]; omega1[n] = omega[n]; 
            // for(int m = 0; m <roundParticle.perimeterPointCount ; m++) 
            // {
            //     xe[n,m] = xp[n] + rp*Mathf.Cos(2.0f*Mathf.PI*m/(float)ne);
            //     ye[n,m] = yp[n] + rp*Mathf.Sin(2.0f*Mathf.PI*m/(float)ne);
            //     ue[n,m] = up[n]  - omega[n]*(ye[n,m] - yp[n]);
            //     ve[n,m] = vp[n]  + omega[n]*(xe[n,m] - xp[n]);
            // } 
            roundParticle.UpdatePerimeter();
        }
        // for(int n = 0; n < particleCount; n++) 
        // { 
        //     fxc[n] = 0.0f;
        //     fyc[n] = 0.0f;
        //     tmp1 = Mathf.Abs(yp[n] + rp); 
        //     if(tmp1 < 2.0f*rp + zeta){
        //         fyc[n] += (yp[n] + rp)*(2.0f*rp - tmp1 + zeta)*(2.0f*rp - tmp1 + zeta)/epsw;
        //     }
        //     tmp1 = Mathf.Abs(xp[n] + rp);
        //     if(tmp1 < 2.0f*rp + zeta){
        //         fxc[n] += (xp[n] + rp)*(2.0f*rp - tmp1 + zeta)*(2.0f*rp - tmp1 + zeta)/epsw;
        //     }
        //     tmp1 = Mathf.Abs((DIM_X-1-xp[n]) + rp);
        //     if(tmp1 < 2.0f*rp + zeta){
        //         fxc[n] -= ((DIM_X-1-xp[n]) + rp)*(2.0f*rp - tmp1 + zeta)*(2.0f*rp - tmp1 + zeta)/epsw;
        //     }

        //     fxp[n] = 0.0f; fyp[n] = 0.0f; torq[n] = 0.0f;
        //     for(int m = 0; m < ne ; m++) 
        //     {
        //         uet[n,m] = 0.0f; vet[n,m] = 0.0f;
        //         // 固体表面の速度を計算
        //         for(int i = (int)xe[n,m] - 3; i < (int)xe[n,m] + 3; i++)
        //         {
        //             for(int j = (int)ye[n,m] - 3; j < (int)ye[n,m] + 3; j++)
        //             {
        //                 tmp1 = Mathf.Abs(xe[n,m] - (float)i);
        //                 tmp2 = Mathf.Abs(ye[n,m] - (float)j);
        //                 if(tmp1 <= 2.0f)
        //                 {
        //                     tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp1/2.0f))/4.0f;
        //                 } 
        //                 else 
        //                 {
        //                     tmp3 = 0.0f;
        //                 }
        //                 if(tmp2 <= 2.0f)
        //                 {
        //                     tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp2/2.0f))/4.0f*tmp3;
        //                 } 
        //                 else 
        //                 {
        //                     tmp3 = 0.0f;
        //                 }
        //                 if((j<DIM_Y&&j>=0) && (i<DIM_X&&i>=0))
        //                 {
        //                     uet[n,m] += u[i,j]*tmp3;
        //                     vet[n,m] += v[i,j]*tmp3;
        //                 }
        //             } 
        //         }
        //         fxe[n,m] = ue[n,m] - uet[n,m];
        //         fye[n,m] = ve[n,m] - vet[n,m];
        //         // 固体が外部に与える力を計算
        //         for(int i = (int)xe[n,m] - 3; i < (int)xe[n,m] + 3; i++)
        //         {
        //             for(int j = (int)ye[n,m] - 3; j < (int)ye[n,m] + 3; j++)
        //             {
        //                 tmp1 = Mathf.Abs(xe[n,m] - (float)i);
        //                 tmp2 = Mathf.Abs(ye[n,m] - (float)j);

        //                 if(tmp1 <= 2.0f)
        //                 {
        //                     tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp1/2.0f))/4.0f;
        //                 } 
        //                 else 
        //                 {
        //                     tmp3 = 0.0f;
        //                 }
        //                 if(tmp2 <= 2.0f)
        //                 {
        //                     tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp2/2.0f))/4.0f*tmp3;
        //                 } 
        //                 else 
        //                 {
        //                     tmp3 = 0.0f;
        //                 }
        //                 if((j<DIM_Y&&j>=0) && (i<DIM_X&&i>=0))
        //                 {
        //                     fx[i,j] += fxe[n,m] * tmp3 * 2.0f*Mathf.PI*rp/(float)ne;
        //                     fy[i,j] += fye[n,m] * tmp3 * 2.0f*Mathf.PI*rp/(float)ne;
        //                 }
        //             } 
        //         }
        //         fxp[n] += fxe[n,m];
        //         fyp[n] += fye[n,m];
        //         torq[n] = torq[n] + fye[n,m]*(xe[n,m] - xp[n]) - fxe[n,m]*(ye[n,m] - yp[n]);
        //     } 
        //     fxp[n] =  -fxp[n]*2.0f*Mathf.PI*rp/(float)ne;
        //     fyp[n] =  -fyp[n]*2.0f*Mathf.PI*rp/(float)ne;
        //     torq[n] = -torq[n]*2.0f*Mathf.PI*rp/(float)ne;
        //     up[n]  = (1.0f + 1.0f/rhop)*up1[n]    - 1.0f/rhop*up2[n]
        //             + (fxp[n] + fxc[n])/rhop/mp + (1.0f - 1.0f/rhop)*gx;
        //     vp[n]  = (1.0f + 1.0f/rhop)*vp1[n]    - 1.0f/rhop*vp2[n]
        //             + (fyp[n] + fyc[n])/rhop/mp + (1.0f - 1.0f/rhop)*gy;
        //     omega[n]= (1.0f + 1.0f/rhop)*omega1[n] - 1.0f/rhop*omega2[n]
        //             + torq[n]/rhop/iip;
        //         xp[n]=    xp[n] + (   up[n] +    up1[n])*0.5f;
        //         yp[n]=    yp[n] + (   vp[n] +    vp1[n])*0.5f;
        //     theta[n]= theta[n] + (omega[n] + omega1[n])*0.5f;
        //         up2[n] =    up1[n];    up1[n] =    up[n];
        //         vp2[n] =    vp1[n];    vp1[n] =    vp[n];
        //     omega2[n] = omega1[n]; omega1[n] = omega[n]; 
        //     for(int m = 0; m <ne ; m++) 
        //     {
        //         xe[n,m] = xp[n] + rp*Mathf.Cos(2.0f*Mathf.PI*m/(float)ne);
        //         ye[n,m] = yp[n] + rp*Mathf.Sin(2.0f*Mathf.PI*m/(float)ne);
        //         ue[n,m] = up[n]  - omega[n]*(ye[n,m] - yp[n]);
        //         ve[n,m] = vp[n]  + omega[n]*(xe[n,m] - xp[n]);
        //     } 
        // }
//         int n,m,i,j;
//             for(n = 0; n < particleCount; n++) { for(m = 0; m <ne ; m++) {
//       uet[n,m] = 0.0f; vet[n,m] = 0.0f;
//     } }

//     for(n = 0; n < particleCount; n++) { for(m = 0; m <ne ; m++) {
//       for(i = (int)xe[n,m] - 3; i < (int)xe[n,m] + 3; i++){
//       for(j = (int)ye[n,m] - 3; j < (int)ye[n,m] + 3; j++){
//         tmp1 = Mathf.Abs(xe[n,m] - (float)i);
//         tmp2 = Mathf.Abs(ye[n,m] - (float)j);

//         if(tmp1 <= 2.0f){
//           tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp1/2.0f))/4.0f;
//         } else {
//           tmp3 = 0.0f;
//         }
//         if(tmp2 <= 2.0f){
//           tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp2/2.0f))/4.0f*tmp3;
//         } else {
//           tmp3 = 0.0f;
//         }
//         if((j<DIM_Y&&j>=0) && (i<DIM_X&&i>=0))
//         {
//             uet[n,m] += u[i,j]*tmp3;
//             vet[n,m] += v[i,j]*tmp3;
//         }
//       } }
//     } }

//     for(n = 0; n < particleCount; n++) { for(m = 0; m <ne ; m++) {
//       fxe[n,m] = 0.0f; fye[n,m] = 0.0f; 
//     } }

//     for(n = 0; n < particleCount; n++) { for(m = 0; m <ne ; m++) {
//       fxe[n,m] = ue[n,m] - uet[n,m];
//       fye[n,m] = ve[n,m] - vet[n,m];
//     } }

// // Spreading: force (fxpe, fype) --> (fx, fy)
//     for(i = 0; i < DIM_X-1; i++) { for(j = 0; j < DIM_Y-1; j++) {
//       fx[i,j] = 0.0f; fy[i,j] = 0.0f;
//     } }

//     for(n = 0; n < particleCount; n++) { for(m = 0; m <ne ; m++) {
//       for(i = (int)xe[n,m] - 3; i < (int)xe[n,m] + 3; i++){
//       for(j = (int)ye[n,m] - 3; j < (int)ye[n,m] + 3; j++){

//         tmp1 = Mathf.Abs(xe[n,m] - (float)i);
//         tmp2 = Mathf.Abs(ye[n,m] - (float)j);

//         if(tmp1 <= 2.0f){
//           tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp1/2.0f))/4.0f;
//         } else {
//           tmp3 = 0.0f;
//         }
//         if(tmp2 <= 2.0f){
//           tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp2/2.0f))/4.0f*tmp3;
//         } else {
//           tmp3 = 0.0f;
//         }
//         if((j<DIM_Y&&j>=0) && (i<DIM_X&&i>=0))
//         {
//             fx[i,j] += fxe[n,m] * tmp3 * 2.0f*Mathf.PI*rp/(float)ne;
//             fy[i,j] += fye[n,m] * tmp3 * 2.0f*Mathf.PI*rp/(float)ne;
//         }
//       } }
//     } }

//     for(n = 0; n < particleCount; n++){
//       fxc[n] = 0.0f;
//       fyc[n] = 0.0f;
//     }

// // particle-wall collision
// // bottom wall 
//     tmp1 = Mathf.Abs(yp[n] + rp); 
//     if(tmp1 < 2.0f*rp + zeta){
//       fyc[n] += (yp[n] + rp)*(2.0f*rp - tmp1 + zeta)*(2.0f*rp - tmp1 + zeta)/epsw;
//     } 

// // particle motion
//     for(n = 0; n < particleCount; n++){
//       fxp[n] = 0.0f; fyp[n] = 0.0f; torq[n] = 0.0f;
//     }
//     for(n = 0; n < particleCount; n++){ for(m = 0; m <ne ; m++){
//        fxp[n] += fxe[n,m];
//        fyp[n] += fye[n,m];
//       torq[n]  = torq[n] + fye[n,m]*(xe[n,m] - xp[n])
//                          - fxe[n,m]*(ye[n,m] - yp[n]);
//     } }
//     for(n = 0; n < particleCount; n++){
//        fxp[n] =  -fxp[n]*2.0f*Mathf.PI*rp/(float)ne;
//        fyp[n] =  -fyp[n]*2.0f*Mathf.PI*rp/(float)ne;
//       torq[n] = -torq[n]*2.0f*Mathf.PI*rp/(float)ne;
//     }

//     for(n = 0; n < particleCount; n++){
//        up[n]  = (1.0f + 1.0f/rhop)*up1[n]    - 1.0f/rhop*up2[n]
//               + (fxp[n] + fxc[n])/rhop/mp + (1.0f - 1.0f/rhop)*gx;

//        vp[n]  = (1.0f + 1.0f/rhop)*vp1[n]    - 1.0f/rhop*vp2[n]
//               + (fyp[n] + fyc[n])/rhop/mp + (1.0f - 1.0f/rhop)*gy;

//       omega[n]= (1.0f + 1.0f/rhop)*omega1[n] - 1.0f/rhop*omega2[n]
//               + torq[n]/rhop/iip;

//          xp[n]=    xp[n] + (   up[n] +    up1[n])*0.5f;
//          yp[n]=    yp[n] + (   vp[n] +    vp1[n])*0.5f;
//       theta[n]= theta[n] + (omega[n] + omega1[n])*0.5f;

//          up2[n] =    up1[n];    up1[n] =    up[n];
//          vp2[n] =    vp1[n];    vp1[n] =    vp[n];
//       omega2[n] = omega1[n]; omega1[n] = omega[n]; 
//     }

//     for(n = 0; n < particleCount; n++) { for(m = 0; m <ne ; m++) {
//       xe[n,m] = xp[n] + rp*Mathf.Cos(2.0f*Mathf.PI*m/(float)ne);
//       ye[n,m] = yp[n] + rp*Mathf.Sin(2.0f*Mathf.PI*m/(float)ne);
//     } } 

//     for(n = 0; n < particleCount; n++) { for(m = 0; m <ne ; m++) {
//       ue[n,m] = up[n]  - omega[n]*(ye[n,m] - yp[n]);
//       ve[n,m] = vp[n]  + omega[n]*(xe[n,m] - xp[n]);
//     } }
    }
}
