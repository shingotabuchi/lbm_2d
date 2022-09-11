using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LBMNatConvComputeShader : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public VectorField vectorField;
    public HeatMapMode mode = HeatMapMode.Speed;
    public InitialTemperatureDistribution initTemp = InitialTemperatureDistribution.XGradient;
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
    // float[] minMaxTempAndSpeed = new float[4];

    float[] cx = new float[9]{0, 1,    0,   -1,    0,     1,    -1,    -1,     1};
    float[] cy = new float[9]{0, 0,    1,    0,   -1,     1,     1,    -1,    -1};
    float[] wg = new float[5]{1f/3f,1f/6f,1f/6f,1f/6f,1f/6f};
    float[] wf = new float[9]{4f/9f,1f/9f,1f/9f,1f/9f,1f/9f,1f/36f,1f/36f,1f/36f,1f/36f};
    float[] rho, u, v, temp, fx, fy,speed;
    float[] f, f0, ftmp;
    float[] g, g0, gtmp;
    
    float umax, umin, tmp, u2, nu, chi, norm, taug, rbetag, h;
    public float pr = 0.71f;
    public float ra =   10000.0f;
    public float tauf = 0.8f;

    public int loopCount = 1;

    public ComputeShader computeShader;
    int kernelInitTempXGradient,kernelCollision,kernelStreaming,kernelAdiabaticBoundaryX,kernelAdiabaticBoundaryY,kernelCalculateTempAndSpeed;
    uint threadGroupSize;
    int threadGroupCount,threadGroupCountDIM_X,threadGroupCountDIM_Y;
    ComputeBuffer rhoComputeBuffer, fComputeBuffer,ftmpComputeBuffer,gComputeBuffer,gtmpComputeBuffer;
    ComputeBuffer uComputeBuffer, vComputeBuffer,tempComputeBuffer,speedComputeBuffer;
    // ComputeBuffer minMaxTempAndSpeedBuffer;
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

        kernelInitTempXGradient = computeShader.FindKernel("InitTempXGradient");
        kernelCollision = computeShader.FindKernel("Collision");
        kernelStreaming = computeShader.FindKernel("Streaming");
        kernelAdiabaticBoundaryX = computeShader.FindKernel("AdiabaticBoundaryX");
        kernelAdiabaticBoundaryY = computeShader.FindKernel("AdiabaticBoundaryY");
        kernelCalculateTempAndSpeed = computeShader.FindKernel("CalculateTempAndSpeed");
        
        computeShader.GetKernelThreadGroupSizes(kernelInitTempXGradient, out threadGroupSize, out _, out _);
        threadGroupCount = (int) ((DIM_X*DIM_Y + (threadGroupSize - 1)) / threadGroupSize);
        threadGroupCountDIM_X = (int) ((DIM_X + (threadGroupSize - 1)) / threadGroupSize);
        threadGroupCountDIM_Y = (int) ((DIM_Y + (threadGroupSize - 1)) / threadGroupSize);
        computeShader.SetInt("DIM_X",DIM_X);
        computeShader.SetInt("DIM_Y",DIM_Y);

        computeShader.SetFloat("tauf",tauf);
        computeShader.SetFloat("taug",taug);
        computeShader.SetFloat("rbetag",rbetag);

        rhoComputeBuffer = new ComputeBuffer(DIM_X*DIM_Y, sizeof(float));
        fComputeBuffer =  new ComputeBuffer(DIM_X*DIM_Y*9, sizeof(float));
        ftmpComputeBuffer =  new ComputeBuffer(DIM_X*DIM_Y*9, sizeof(float));
        gComputeBuffer =  new ComputeBuffer(DIM_X*DIM_Y*5, sizeof(float));
        gtmpComputeBuffer =  new ComputeBuffer(DIM_X*DIM_Y*5, sizeof(float));
        uComputeBuffer =  new ComputeBuffer(DIM_X*DIM_Y, sizeof(float));
        vComputeBuffer =  new ComputeBuffer(DIM_X*DIM_Y, sizeof(float));
        tempComputeBuffer = new ComputeBuffer(DIM_X*DIM_Y, sizeof(float));
        speedComputeBuffer = new ComputeBuffer(DIM_X*DIM_Y, sizeof(float));
        // minMaxTempAndSpeedBuffer = new ComputeBuffer(4, sizeof(float));

        computeShader.SetBuffer(kernelInitTempXGradient, "rho", rhoComputeBuffer);
        computeShader.SetBuffer(kernelInitTempXGradient, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelInitTempXGradient, "g", gComputeBuffer);
        computeShader.SetBuffer(kernelInitTempXGradient, "u", uComputeBuffer);
        computeShader.SetBuffer(kernelInitTempXGradient, "v", vComputeBuffer);
        computeShader.SetBuffer(kernelInitTempXGradient, "temp", tempComputeBuffer);
        computeShader.SetBuffer(kernelInitTempXGradient, "speed", speedComputeBuffer);

        computeShader.SetBuffer(kernelCollision, "rho", rhoComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "ftmp", ftmpComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "g", gComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "gtmp", gtmpComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "u", uComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "v", vComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "temp", tempComputeBuffer);

        computeShader.SetBuffer(kernelStreaming, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelStreaming, "ftmp", ftmpComputeBuffer);
        computeShader.SetBuffer(kernelStreaming, "g", gComputeBuffer);
        computeShader.SetBuffer(kernelStreaming, "gtmp", gtmpComputeBuffer);

        computeShader.SetBuffer(kernelAdiabaticBoundaryX, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelAdiabaticBoundaryX, "g", gComputeBuffer);
        computeShader.SetBuffer(kernelAdiabaticBoundaryY, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelAdiabaticBoundaryY, "g", gComputeBuffer);

        computeShader.SetBuffer(kernelCalculateTempAndSpeed, "u", uComputeBuffer);
        computeShader.SetBuffer(kernelCalculateTempAndSpeed, "v", vComputeBuffer);
        computeShader.SetBuffer(kernelCalculateTempAndSpeed, "rho", rhoComputeBuffer);
        computeShader.SetBuffer(kernelCalculateTempAndSpeed, "temp", tempComputeBuffer);
        computeShader.SetBuffer(kernelCalculateTempAndSpeed, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelCalculateTempAndSpeed, "g", gComputeBuffer);
        computeShader.SetBuffer(kernelCalculateTempAndSpeed, "speed", speedComputeBuffer);
        // computeShader.SetBuffer(kernelCalculateTempAndSpeed, "minMaxTempAndSpeed", minMaxTempAndSpeedBuffer);

        computeShader.Dispatch(kernelInitTempXGradient, threadGroupCount, 1, 1);

        speed = new float[DIM_X*DIM_Y];
        u = new float[DIM_X*DIM_Y];
        v = new float[DIM_X*DIM_Y];
        temp = new float[DIM_X*DIM_Y];
    }

    void LBMStep()
    {
        Collision();
        // Force();
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
        // minMaxTempAndSpeedBuffer.GetData(minMaxTempAndSpeed);
        // minTemp = minMaxTempAndSpeed[0];
        // maxTemp = minMaxTempAndSpeed[1];
        // minSpeed = minMaxTempAndSpeed[2];
        // maxSpeed = minMaxTempAndSpeed[3];
        speedComputeBuffer.GetData(speed);
        tempComputeBuffer.GetData(temp);
        uComputeBuffer.GetData(u);
        vComputeBuffer.GetData(v);
        maxSpeed = 0f;
        minSpeed = Mathf.Infinity;
        maxTemp = 0f;
        minTemp = Mathf.Infinity;
        for (int i = 0; i < plotPixels.Length; i++)
        {
            maxSpeed = Mathf.Max(speed[i%DIM_X+(i/DIM_X)*DIM_X],maxSpeed);
            minSpeed = Mathf.Min(speed[i%DIM_X+(i/DIM_X)*DIM_X],minSpeed);
            maxTemp = Mathf.Max(temp[i%DIM_X+(i/DIM_X)*DIM_X],maxTemp);
            minTemp = Mathf.Min(temp[i%DIM_X+(i/DIM_X)*DIM_X],minTemp);
        }
        for (int i = 0; i < plotPixels.Length; i++)
        {
            if(normalizeHeatMap)
            {
                if(mode == HeatMapMode.Speed) 
                plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X+(i/DIM_X)*DIM_X]-minSpeed,maxSpeed-minSpeed);
                else
                plotPixels[i] = colorHeatMap.GetColorForValue(temp[i%DIM_X+(i/DIM_X)*DIM_X]-minTemp,maxTemp-minTemp);
            }
            else
            {
                if(mode == HeatMapMode.Speed) 
                plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X+(i/DIM_X)*DIM_X],maxSpeed);
                else
                plotPixels[i] = colorHeatMap.GetColorForValue(temp[i%DIM_X+(i/DIM_X)*DIM_X],maxTemp);
            }
        }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
        vectorField.UpdateVectors(u,v,DIM_Y,maxSpeed);
    }

    void Collision()
    {
        // computeShader.SetBuffer(kernelCollision, "rho", rhoComputeBuffer);
        // computeShader.SetBuffer(kernelCollision, "f", fComputeBuffer);
        // computeShader.SetBuffer(kernelCollision, "ftmp", ftmpComputeBuffer);
        // computeShader.SetBuffer(kernelCollision, "g", gComputeBuffer);
        // computeShader.SetBuffer(kernelCollision, "gtmp", gtmpComputeBuffer);
        // computeShader.SetBuffer(kernelCollision, "u", uComputeBuffer);
        // computeShader.SetBuffer(kernelCollision, "v", vComputeBuffer);
        // computeShader.SetBuffer(kernelCollision, "temp", tempComputeBuffer);
        
        computeShader.Dispatch(kernelCollision, threadGroupCount, 1, 1);
        

    }

    // void Force()
    // {
    //     for (int i = 0; i < DIM_X; i++)
    //     {
    //         for (int j = 0; j < DIM_Y; j++)
    //         {
    //             fx[i + j*DIM_X] = 0.0f; fy[i + j*DIM_X] = rbetag*(temp[i + j*DIM_X] - 0.5f);
    //             for (int k = 0; k < 9; k++)
    //             {
    //                 f[k + (i + j*DIM_X)*9] = f[k + (i + j*DIM_X)*9] + 3f*wf[k]*(cx[k]*fx[i + j*DIM_X] + cy[k]*fy[i + j*DIM_X]);
    //                 ftmp[k + (i + j*DIM_X)*9] = f[k + (i + j*DIM_X)*9];
    //             }
    //         }
    //     }
    // }
    
    void Streaming()
    {
        // ftmp = (float[,,])(f.Clone());
        // gtmp = (float[,,])(g.Clone());

        // for(int i = 0; i < DIM_X; i++)
        // { 
        //     for(int j = 0; j < DIM_Y; j++)
        //     { 
        //         for(int k = 0; k < 9; k++)
        //         {
        //             int im = i + (int)cx[k]; 
        //             int jm = j + (int)cy[k];
        //             if((jm!=DIM_Y&&jm!=-1) && (im!=DIM_X&&im!=-1))
        //             {
        //                 f[k + (im + jm*DIM_X)*9] = ftmp[k + (i + j*DIM_X)*9];
        //                 if(k<5)
        //                 g[k + (im + jm*DIM_X)*5] = gtmp[k + (i + j*DIM_X)*5];
        //             }
        //         } 
        //     }
        // }

        // computeShader.SetBuffer(kernelStreaming, "f", fComputeBuffer);
        // computeShader.SetBuffer(kernelStreaming, "ftmp", ftmpComputeBuffer);
        // computeShader.SetBuffer(kernelStreaming, "g", gComputeBuffer);
        // computeShader.SetBuffer(kernelStreaming, "gtmp", gtmpComputeBuffer);
        computeShader.Dispatch(kernelStreaming, threadGroupCount, 1, 1);
        // fComputeBuffer.GetData(f);
        // gComputeBuffer.GetData(g);
    }

    void Boundaries()
    {
        // computeShader.SetBuffer(kernelAdiabaticBoundaryX, "f", fComputeBuffer);
        // computeShader.SetBuffer(kernelAdiabaticBoundaryX, "g", gComputeBuffer);
        // computeShader.SetBuffer(kernelAdiabaticBoundaryY, "f", fComputeBuffer);
        // computeShader.SetBuffer(kernelAdiabaticBoundaryY, "g", gComputeBuffer);
        
        computeShader.Dispatch(kernelAdiabaticBoundaryX, threadGroupCountDIM_X, 1, 1);
        computeShader.Dispatch(kernelAdiabaticBoundaryY, threadGroupCountDIM_Y, 1, 1);
        // fComputeBuffer.GetData(f);
        // gComputeBuffer.GetData(g);
    }

    void UpdateSpeedAndTemperature()
    {
        maxSpeed = 0f;
        minSpeed = Mathf.Infinity;
        maxTemp = 0f;
        minTemp = Mathf.Infinity;
        // minMaxTempAndSpeed[0] = minTemp;
        // minMaxTempAndSpeed[1] = maxTemp;
        // minMaxTempAndSpeed[2] = minSpeed;
        // minMaxTempAndSpeed[3] = maxSpeed;
        // minMaxTempAndSpeedBuffer.SetData(minMaxTempAndSpeed);
        // computeShader.SetBuffer(kernelCalculateTempAndSpeed,"minMaxTempAndSpeed",minMaxTempAndSpeedBuffer);
        computeShader.Dispatch(kernelCalculateTempAndSpeed,threadGroupCount,1,1);
        // for(int i = 0; i < DIM_X; i++)
        // { 
        //     for(int j = 0; j < DIM_Y; j++)
        //     {
        //         u[i + j*DIM_X] = 0f; v[i + j*DIM_X] = 0f;
        //         rho[i + j*DIM_X] = f[0+(i+j*DIM_X)*9]; 
        //         temp[i + j*DIM_X] =  g[0+(i+j*DIM_X)*5];
        //         for(int k = 1; k <= 8; k++)
        //         {
        //             rho[i + j*DIM_X] = rho[i + j*DIM_X] + f[k + (i + j*DIM_X)*9];
        //             u[i + j*DIM_X] =   u[i + j*DIM_X] + f[k + (i + j*DIM_X)*9]*cx[k];
        //             v[i + j*DIM_X] =   v[i + j*DIM_X] + f[k + (i + j*DIM_X)*9]*cy[k];
        //             if(k<5){
        //             temp[i + j*DIM_X] =   temp[i + j*DIM_X] + g[k + (i + j*DIM_X)*5];}
                    
        //         } 
        //         u[i + j*DIM_X] = u[i + j*DIM_X]/rho[i + j*DIM_X];
        //         v[i + j*DIM_X] = v[i + j*DIM_X]/rho[i + j*DIM_X];
        //         speed[i + j*DIM_X] = Mathf.Sqrt(u[i + j*DIM_X]*u[i + j*DIM_X] + v[i + j*DIM_X]*v[i + j*DIM_X]);
                
        //         maxTemp = Mathf.Max(maxTemp,temp[i + j*DIM_X]);
        //         minTemp = Mathf.Min(minTemp,temp[i + j*DIM_X]);
        //         maxSpeed = Mathf.Max(maxSpeed,speed[i + j*DIM_X]);
        //         minSpeed = Mathf.Min(minSpeed,speed[i + j*DIM_X]);
        //         maxRho = Mathf.Max(maxRho,rho[i + j*DIM_X]);
        //         minRho = Mathf.Min(minRho,rho[i + j*DIM_X]);
        //     } 
        // }
        // rhoComputeBuffer.SetData(rho);
        // fComputeBuffer.SetData(f);
        // gComputeBuffer.SetData(g);
        // uComputeBuffer.SetData(u);
        // vComputeBuffer.SetData(v);
        // tempComputeBuffer.SetData(temp);
    }
}