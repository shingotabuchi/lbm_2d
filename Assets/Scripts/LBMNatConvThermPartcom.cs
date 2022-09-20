using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LBMNatConvThermPartcom : MonoBehaviour
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

    float setWallTemp1,setWallTemp2,setWallTemp3,setWallTemp4;
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

    float[] rho, u, v, temp,speed,fx,fy;
    float[] f, f0, ftmp;
    float[] g, g0, gtmp;

    // float[,] rho, u, v, temp, forceFromGravityX, forceFromGravityY,speed, fx,fy;
    // float[,,] f, f0, ftmp;
    // float[,,] g, g0, gtmp;

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
    RoundParticlesForCompute roundParticles;
    
    public InitialTemperatureDistribution initTemp = InitialTemperatureDistribution.XGradient;

    public GameObject particlePrefab;
    public Transform particleParent;
    Vector2 plotSizeDelta; 
    float plotScale; 
    [Range(0.0f, 1.0f)]
    public float particleTemp = 1f;
    [Range(0.0f, 1.0f)]
    public float particleLowerTemp = 0f;

    public ComputeShader computeShader;
    int kernelInitTempXGradient,kernelInitTempYGradient,kernelInitTempAllZero;
    int kernelCollision,kernelStreaming;
    int kernelAdiabaticBoundaryX,
        kernelAdiabaticBoundaryY,
        kernelConstantBoundaryX,
        kernelConstantBoundaryY,
        kernelCalculateTempAndSpeed;
    int kernelConstantAdiabaticBoundaryX,
        kernelConstantAdiabaticBoundaryY,
        kernelAdiabaticConstantBoundaryX,
        kernelAdiabaticConstantBoundaryY;
    uint threadGroupSize;
    int threadGroupCount,threadGroupCountDIM_X,threadGroupCountDIM_Y;
    ComputeBuffer rhoComputeBuffer, fComputeBuffer,ftmpComputeBuffer,gComputeBuffer,gtmpComputeBuffer;
    ComputeBuffer uComputeBuffer, vComputeBuffer,tempComputeBuffer,speedComputeBuffer;
    ComputeBuffer fxComputeBuffer, fyComputeBuffer;
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
        gravity[0] = 0f;
        gravity[1] = -rbetag/beta;
        
        float[] spawnPos = new float[particleCount*2];
        GameObject[] objs = new GameObject[particleCount];
        float scale = (2f*particleRadius*plotImage.transform.GetComponent<RectTransform>().sizeDelta.x * plotImage.transform.localScale.x)/(DIM_X*particlePrefab.GetComponent<RectTransform>().sizeDelta.x);
        plotSizeDelta = plotImage.transform.GetComponent<RectTransform>().sizeDelta;
        plotScale = plotImage.transform.localScale.x;
        for (int i = 0; i < particleCount; i++)
        {
            spawnPos[i + 0*particleCount] = UnityEngine.Random.Range(particleRadius,DIM_X-particleRadius);
            spawnPos[i + 1*particleCount] = UnityEngine.Random.Range(particleRadius,DIM_Y-particleRadius);
            objs[i] = Instantiate(particlePrefab,particleParent);
            objs[i].GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
                                                            + new Vector2((spawnPos[i + 0*particleCount]*plotSizeDelta.x)/DIM_X,(spawnPos[i + 1*particleCount]*plotSizeDelta.y)/DIM_Y))*plotScale;
            objs[i].transform.localScale = new Vector3(scale,scale,1);
        }
        roundParticles = new RoundParticlesForCompute(particleCount,particleDensity,particleRadius,spawnPos,objs,Random.Range(0f,2f*Mathf.PI));

        kernelInitTempXGradient = computeShader.FindKernel("InitTempXGradient");
        kernelInitTempYGradient = computeShader.FindKernel("InitTempYGradient");
        kernelInitTempAllZero = computeShader.FindKernel("InitTempAllZero");
        kernelCollision = computeShader.FindKernel("Collision");
        kernelStreaming = computeShader.FindKernel("Streaming");
        kernelAdiabaticBoundaryX = computeShader.FindKernel("AdiabaticBoundaryX");
        kernelAdiabaticBoundaryY = computeShader.FindKernel("AdiabaticBoundaryY");
        kernelConstantBoundaryX = computeShader.FindKernel("ConstantBoundaryX");
        kernelConstantBoundaryY = computeShader.FindKernel("ConstantBoundaryY");
        kernelConstantAdiabaticBoundaryX = computeShader.FindKernel("ConstantAdiabaticBoundaryX");
        kernelConstantAdiabaticBoundaryY = computeShader.FindKernel("ConstantAdiabaticBoundaryY");
        kernelAdiabaticConstantBoundaryX = computeShader.FindKernel("AdiabaticConstantBoundaryX");
        kernelAdiabaticConstantBoundaryY = computeShader.FindKernel("AdiabaticConstantBoundaryY");
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

        computeShader.SetFloat("wallTemp1",wallTemp1);
        computeShader.SetFloat("wallTemp2",wallTemp2);
        computeShader.SetFloat("wallTemp3",wallTemp3);
        computeShader.SetFloat("wallTemp4",wallTemp4);

        setWallTemp1 = wallTemp1;
        setWallTemp2 = wallTemp2;
        setWallTemp3 = wallTemp3;
        setWallTemp4 = wallTemp4;

        rhoComputeBuffer = new ComputeBuffer(DIM_X*DIM_Y, sizeof(float));
        fComputeBuffer =  new ComputeBuffer(DIM_X*DIM_Y*9, sizeof(float));
        ftmpComputeBuffer =  new ComputeBuffer(DIM_X*DIM_Y*9, sizeof(float));
        gComputeBuffer =  new ComputeBuffer(DIM_X*DIM_Y*5, sizeof(float));
        gtmpComputeBuffer =  new ComputeBuffer(DIM_X*DIM_Y*5, sizeof(float));
        uComputeBuffer =  new ComputeBuffer(DIM_X*DIM_Y, sizeof(float));
        vComputeBuffer =  new ComputeBuffer(DIM_X*DIM_Y, sizeof(float));
        fxComputeBuffer =  new ComputeBuffer(DIM_X*DIM_Y, sizeof(float));
        fyComputeBuffer =  new ComputeBuffer(DIM_X*DIM_Y, sizeof(float));
        tempComputeBuffer = new ComputeBuffer(DIM_X*DIM_Y, sizeof(float));
        speedComputeBuffer = new ComputeBuffer(DIM_X*DIM_Y, sizeof(float));
        // minMaxTempAndSpeedBuffer = new ComputeBuffer(4, sizeof(float));

        if(initTemp == InitialTemperatureDistribution.XGradient)
        {
            computeShader.SetBuffer(kernelInitTempXGradient, "rho", rhoComputeBuffer);
            computeShader.SetBuffer(kernelInitTempXGradient, "f", fComputeBuffer);
            computeShader.SetBuffer(kernelInitTempXGradient, "g", gComputeBuffer);
            computeShader.SetBuffer(kernelInitTempXGradient, "u", uComputeBuffer);
            computeShader.SetBuffer(kernelInitTempXGradient, "v", vComputeBuffer);
            computeShader.SetBuffer(kernelInitTempXGradient, "px", fxComputeBuffer);
            computeShader.SetBuffer(kernelInitTempXGradient, "py", fyComputeBuffer);
            computeShader.SetBuffer(kernelInitTempXGradient, "temp", tempComputeBuffer);
            computeShader.SetBuffer(kernelInitTempXGradient, "speed", speedComputeBuffer);
            computeShader.Dispatch(kernelInitTempXGradient, threadGroupCount, 1, 1);
        }
        if(initTemp == InitialTemperatureDistribution.YGradient)
        {
            computeShader.SetBuffer(kernelInitTempYGradient, "rho", rhoComputeBuffer);
            computeShader.SetBuffer(kernelInitTempYGradient, "f", fComputeBuffer);
            computeShader.SetBuffer(kernelInitTempYGradient, "g", gComputeBuffer);
            computeShader.SetBuffer(kernelInitTempYGradient, "u", uComputeBuffer);
            computeShader.SetBuffer(kernelInitTempYGradient, "v", vComputeBuffer);
            computeShader.SetBuffer(kernelInitTempYGradient, "temp", tempComputeBuffer);
            computeShader.SetBuffer(kernelInitTempYGradient, "speed", speedComputeBuffer);
            computeShader.SetBuffer(kernelInitTempYGradient, "px", fxComputeBuffer);
            computeShader.SetBuffer(kernelInitTempYGradient, "py", fyComputeBuffer);
            computeShader.Dispatch(kernelInitTempYGradient, threadGroupCount, 1, 1);
        }
        if(initTemp == InitialTemperatureDistribution.AllZero)
        {
            computeShader.SetBuffer(kernelInitTempAllZero, "rho", rhoComputeBuffer);
            computeShader.SetBuffer(kernelInitTempAllZero, "f", fComputeBuffer);
            computeShader.SetBuffer(kernelInitTempAllZero, "g", gComputeBuffer);
            computeShader.SetBuffer(kernelInitTempAllZero, "u", uComputeBuffer);
            computeShader.SetBuffer(kernelInitTempAllZero, "v", vComputeBuffer);
            computeShader.SetBuffer(kernelInitTempAllZero, "temp", tempComputeBuffer);
            computeShader.SetBuffer(kernelInitTempAllZero, "speed", speedComputeBuffer);
            computeShader.SetBuffer(kernelInitTempAllZero, "px", fxComputeBuffer);
            computeShader.SetBuffer(kernelInitTempAllZero, "py", fyComputeBuffer);
            computeShader.Dispatch(kernelInitTempAllZero, threadGroupCount, 1, 1);
        }
        
        computeShader.SetBuffer(kernelCollision, "rho", rhoComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "ftmp", ftmpComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "g", gComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "gtmp", gtmpComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "u", uComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "v", vComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "temp", tempComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "px", fxComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "py", fyComputeBuffer);

        computeShader.SetBuffer(kernelStreaming, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelStreaming, "ftmp", ftmpComputeBuffer);
        computeShader.SetBuffer(kernelStreaming, "g", gComputeBuffer);
        computeShader.SetBuffer(kernelStreaming, "gtmp", gtmpComputeBuffer);

        computeShader.SetBuffer(kernelAdiabaticBoundaryX, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelAdiabaticBoundaryX, "g", gComputeBuffer);
        computeShader.SetBuffer(kernelAdiabaticBoundaryY, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelAdiabaticBoundaryY, "g", gComputeBuffer);
        computeShader.SetBuffer(kernelConstantBoundaryX, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelConstantBoundaryX, "g", gComputeBuffer);
        computeShader.SetBuffer(kernelConstantBoundaryY, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelConstantBoundaryY, "g", gComputeBuffer);

        computeShader.SetBuffer(kernelConstantAdiabaticBoundaryX, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelConstantAdiabaticBoundaryX, "g", gComputeBuffer);
        computeShader.SetBuffer(kernelConstantAdiabaticBoundaryY, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelConstantAdiabaticBoundaryY, "g", gComputeBuffer);
        computeShader.SetBuffer(kernelAdiabaticConstantBoundaryX, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelAdiabaticConstantBoundaryX, "g", gComputeBuffer);
        computeShader.SetBuffer(kernelAdiabaticConstantBoundaryY, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelAdiabaticConstantBoundaryY, "g", gComputeBuffer);

        computeShader.SetBuffer(kernelCalculateTempAndSpeed, "u", uComputeBuffer);
        computeShader.SetBuffer(kernelCalculateTempAndSpeed, "v", vComputeBuffer);
        computeShader.SetBuffer(kernelCalculateTempAndSpeed, "rho", rhoComputeBuffer);
        computeShader.SetBuffer(kernelCalculateTempAndSpeed, "temp", tempComputeBuffer);
        computeShader.SetBuffer(kernelCalculateTempAndSpeed, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelCalculateTempAndSpeed, "g", gComputeBuffer);
        computeShader.SetBuffer(kernelCalculateTempAndSpeed, "speed", speedComputeBuffer);

        speed = new float[DIM_X*DIM_Y];
        u = new float[DIM_X*DIM_Y];
        v = new float[DIM_X*DIM_Y];
        fx = new float[DIM_X*DIM_Y];
        fy = new float[DIM_X*DIM_Y];
        for (int i = 0; i < DIM_X; i++)
        {
            for (int j = 0; j < DIM_Y; j++)
            {
                fx[i + j*DIM_X] = 0f;
                fy[i + j*DIM_X] = 0f;
            }
        }
        temp = new float[DIM_X*DIM_Y];
    }

    void LBMStep()
    {
        // fxComputeBuffer.SetData(fx);
        // fyComputeBuffer.SetData(fy);
        computeShader.Dispatch(kernelCollision, threadGroupCount, 1, 1);
        computeShader.Dispatch(kernelStreaming, threadGroupCount, 1, 1);
        // computeShader.Dispatch(kernelAdiabaticBoundaryX, threadGroupCountDIM_X, 1, 1);
        // computeShader.Dispatch(kernelAdiabaticBoundaryY, threadGroupCountDIM_Y, 1, 1);
        if(wallboundaries[2] == BoundaryType.Adiabatic && wallboundaries[3] == BoundaryType.Adiabatic) 
        computeShader.Dispatch(kernelAdiabaticBoundaryX, threadGroupCountDIM_X, 1, 1);
        if(wallboundaries[2] == BoundaryType.Constant && wallboundaries[3] == BoundaryType.Constant) 
        computeShader.Dispatch(kernelConstantBoundaryX, threadGroupCountDIM_X, 1, 1);
        if(wallboundaries[0] == BoundaryType.Adiabatic && wallboundaries[1] == BoundaryType.Adiabatic) 
        computeShader.Dispatch(kernelAdiabaticBoundaryY, threadGroupCountDIM_Y, 1, 1);
        if(wallboundaries[0] == BoundaryType.Constant && wallboundaries[1] == BoundaryType.Constant)  
        computeShader.Dispatch(kernelConstantBoundaryY, threadGroupCountDIM_Y, 1, 1);

        if(wallboundaries[2] == BoundaryType.Constant && wallboundaries[3] == BoundaryType.Adiabatic) 
        computeShader.Dispatch(kernelConstantAdiabaticBoundaryX, threadGroupCountDIM_X, 1, 1);
        if(wallboundaries[2] == BoundaryType.Adiabatic && wallboundaries[3] == BoundaryType.Constant) 
        computeShader.Dispatch(kernelAdiabaticConstantBoundaryX, threadGroupCountDIM_X, 1, 1);
        if(wallboundaries[0] == BoundaryType.Constant && wallboundaries[1] == BoundaryType.Adiabatic) 
        computeShader.Dispatch(kernelConstantAdiabaticBoundaryY, threadGroupCountDIM_Y, 1, 1);
        if(wallboundaries[0] == BoundaryType.Adiabatic && wallboundaries[1] == BoundaryType.Constant)  
        computeShader.Dispatch(kernelAdiabaticConstantBoundaryY, threadGroupCountDIM_Y, 1, 1);


        computeShader.Dispatch(kernelCalculateTempAndSpeed,threadGroupCount,1,1);
        uComputeBuffer.GetData(u);
        vComputeBuffer.GetData(v);
        // fxComputeBuffer.GetData(fx);
        // fyComputeBuffer.GetData(fy);
        // ImmersedBoundary();
    }

    // Update is called once per frame
    void Update()
    {
        if(setWallTemp1!=wallTemp1){computeShader.SetFloat("wallTemp1",wallTemp1);setWallTemp1 = wallTemp1;}
        if(setWallTemp2!=wallTemp2){computeShader.SetFloat("wallTemp2",wallTemp2);setWallTemp2 = wallTemp2;}
        if(setWallTemp3!=wallTemp3){computeShader.SetFloat("wallTemp3",wallTemp3);setWallTemp3 = wallTemp3;}
        if(setWallTemp4!=wallTemp4){computeShader.SetFloat("wallTemp4",wallTemp4);setWallTemp4 = wallTemp4;}
        for (int i = 0; i < loopCount; i++)
        {
            LBMStep();
        }
        UpdatePlot();
    }

    void UpdatePlot()
    {
        speedComputeBuffer.GetData(speed);
        tempComputeBuffer.GetData(temp);
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
        if(maxSpeed == 0f)maxSpeed = 1f;
        if(maxTemp == 0f)maxTemp = 1f;
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
        // for (int i = 0; i < particleCount; i++)
        // {
        //     roundParticles[i].obj.GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
        //                                                                         + new Vector2((roundParticles[i].pos[0]*plotSizeDelta.x)/DIM_X,(roundParticles[i].pos[1]*plotSizeDelta.y)/DIM_Y))*plotScale;
        //     roundParticles[i].obj.transform.rotation = Quaternion.Euler(0,0,((roundParticles[i].theta*180f)/Mathf.PI));
        // }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
        vectorField.UpdateVectors(u,v,DIM_Y,maxSpeed);
    }

    void ImmersedBoundary()
    {
        float tmp1,tmp2,tmp3;
        for(int n = 0; n < particleCount; n++) 
        { 
            roundParticles.forceFromCollisions[n + 0*particleCount] = 0f;
            roundParticles.forceFromCollisions[n + 1*particleCount] = 0f;

            for (int k = 0; k < particleCount; k++)
            {
                if(k==n) continue;
                for (int i = 0; i < 2; i++)
                {
                    tmp1 = roundParticles.ParticleDistance(n,k);
                    if(tmp1 < 2.0f*roundParticles.radius[n] + zeta)
                    {
                        roundParticles.forceFromCollisions[n + i*particleCount] += (roundParticles.pos[n + i*particleCount] - roundParticles.pos[k + i*particleCount])*(2.0f*roundParticles.radius[n] - tmp1 + zeta)*(2.0f*roundParticles.radius[n] - tmp1 + zeta)/epsw;
                    }
                }
            }

            float wallDistance;
            wallDistance = Mathf.Abs(roundParticles.pos[n + 1*particleCount]);
            if(wallDistance < roundParticles.radius[n]) roundParticles.pos[n + 1*particleCount] = roundParticles.radius[n];
            wallDistance = Mathf.Abs(DIM_Y-1-roundParticles.pos[n + 1*particleCount]);
            if(wallDistance < roundParticles.radius[n]) roundParticles.pos[n + 1*particleCount] = DIM_Y-1-roundParticles.radius[n];
            wallDistance = Mathf.Abs(roundParticles.pos[n + 0*particleCount]);
            if(wallDistance < roundParticles.radius[n]) roundParticles.pos[n + 0*particleCount] =roundParticles.radius[n];
            wallDistance = Mathf.Abs(DIM_X-1-roundParticles.pos[n + 0*particleCount]);
            if(wallDistance < roundParticles.radius[n]) roundParticles.pos[n + 0*particleCount] = DIM_X-1-roundParticles.radius[n];

            roundParticles.forceFromFluid[n + 0*particleCount] = 0f;
            roundParticles.forceFromFluid[n + 1*particleCount] = 0f;
            roundParticles.torque[n] = 0f;
            for(int m = 0; m < roundParticles.perimeterPointCount[n] ; m++) 
            {
                roundParticles.perimeterFluidVel[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] = 0f;
                roundParticles.perimeterFluidVel[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] = 0f;
                float angle = Vector2.Angle(
                    new Vector2(roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] - roundParticles.pos[n + 0*particleCount],roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] - roundParticles.pos[n + 1*particleCount]),
                    new Vector2(-Mathf.Sin(roundParticles.theta[n]),Mathf.Cos(roundParticles.theta[n]))
                );
                if((int)roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] >= 0 && (int)roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] < DIM_X
                    && (int)roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] >= 0 && (int)roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] < DIM_Y
                )
                {
                    if(angle>90f)temp[(int)roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] + ((int)roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount])*DIM_X] = particleTemp;
                    // else temp[(int)roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount],(int)roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount]] = particleLowerTemp;
                }
                // else temp[(int)roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount],(int)roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount]] = particleLowerTemp;
                // temp[(int)roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount],(int)roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount]] = 1f;
                // 固体表面の速度を計算
                for(int i = (int)roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] - 3; i < (int)roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] + 3; i++)
                {
                    for(int j = (int)roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] - 3; j < (int)roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] + 3; j++)
                    {
                        tmp1 = Mathf.Abs(roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] - (float)i);
                        tmp2 = Mathf.Abs(roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] - (float)j);
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
                            roundParticles.perimeterFluidVel[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] += u[i + j*DIM_X]*tmp3;
                            roundParticles.perimeterFluidVel[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] += v[i + j*DIM_X]*tmp3;
                        }
                    } 
                }
                roundParticles.forceOnPerimeter[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] = roundParticles.perimeterVel[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] - roundParticles.perimeterFluidVel[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount];
                roundParticles.forceOnPerimeter[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] = roundParticles.perimeterVel[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] - roundParticles.perimeterFluidVel[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount];

                // 固体が外部に与える力を計算
                for(int i = (int)roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] - 3; i < (int)roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] + 3; i++)
                {
                    for(int j = (int)roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] - 3; j < (int)roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] + 3; j++)
                    {
                        tmp1 = Mathf.Abs(roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] - (float)i);
                        tmp2 = Mathf.Abs(roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] - (float)j);
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
                            fx[i + j*DIM_X] += roundParticles.forceOnPerimeter[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] * tmp3 * 2.0f*Mathf.PI*roundParticles.radius[n]/(float)roundParticles.perimeterPointCount[n];
                            fy[i + j*DIM_X] += roundParticles.forceOnPerimeter[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] * tmp3 * 2.0f*Mathf.PI*roundParticles.radius[n]/(float)roundParticles.perimeterPointCount[n];
                        }
                    } 
                }
                roundParticles.forceFromFluid[n + 0*particleCount] += roundParticles.forceOnPerimeter[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount];
                roundParticles.forceFromFluid[n + 1*particleCount] += roundParticles.forceOnPerimeter[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount];
                roundParticles.torque[n] += roundParticles.forceOnPerimeter[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] * (roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] - roundParticles.pos[n + 0*particleCount]) 
                                        - roundParticles.forceOnPerimeter[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] * (roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] - roundParticles.pos[n + 1*particleCount]);
            } 

            roundParticles.forceFromFluid[n + 0*particleCount] *= -2f*Mathf.PI*roundParticles.radius[n]/(float)roundParticles.perimeterPointCount[n];  
            roundParticles.forceFromFluid[n + 1*particleCount] *= -2f*Mathf.PI*roundParticles.radius[n]/(float)roundParticles.perimeterPointCount[n];  
            roundParticles.torque[n] *= -2f*Mathf.PI*roundParticles.radius[n]/(float)roundParticles.perimeterPointCount[n];  

            roundParticles.UpdatePosVel(gravity);
            roundParticles.UpdateOmegaTheta();
            roundParticles.UpdatePerimeter();
        }
    }
}
