using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class D2Q4ComputeShader : MonoBehaviour
{
    public Image plotImage;
    public int resInt;
    [Range(0.25f, 2.0f)]
    public float alpha;
    public BoundaryType[] wallboundaries = new BoundaryType[]{BoundaryType.Constant,BoundaryType.Constant,BoundaryType.Constant,BoundaryType.Constant};
    [Range(0.0f, 1.0f)]
    public float wallTemp1 = 1f;
    [Range(0.0f, 1.0f)]
    public float wallTemp2 = 0f;
    [Range(0.0f, 1.0f)]
    public float wallTemp3 = 0f;
    [Range(0.0f, 1.0f)]
    public float wallTemp4 = 0f;

    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    // float[] f1,f2,f3,f4;
    float[] f,rho,ftmp;
    float dx = 1f;
    float dt = 1f;
    float csq,omega,feq;
    
    public int loopCount = 1;

    public ComputeShader computeShader;
    ComputeBuffer rhoComputeBuffer, fComputeBuffer;
    // ComputeBuffer f1ComputeBuffer,f2ComputeBuffer,f3ComputeBuffer,f4ComputeBuffer;
    ComputeBuffer ftmpComputeBuffer;
    // ComputeBuffer f1tmpComputeBuffer;
    // ComputeBuffer f2tmpComputeBuffer;
    // ComputeBuffer f3tmpComputeBuffer;
    // ComputeBuffer f4tmpComputeBuffer;
    int kernelInit,kernelCollision,kernelStreaming,kernelConstantBoundary,kernelGetRho;
    uint threadGroupSize;
    int threadGroups;
    int threadGroupsResInt;

    void Start()
    {  
        plotTexture = new Texture2D(resInt,resInt);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,resInt,resInt),Vector2.zero);
        csq=dx*dx/(dt*dt);

        kernelInit = computeShader.FindKernel("D2Q4Init");
        kernelCollision = computeShader.FindKernel("D2Q4Collision");
        kernelStreaming = computeShader.FindKernel("D2Q4Streaming");
        kernelConstantBoundary = computeShader.FindKernel("D2Q4ConstantBoundary");
        kernelGetRho = computeShader.FindKernel("D2Q4GetRho");

        computeShader.GetKernelThreadGroupSizes(kernelInit, out threadGroupSize, out _, out _);
        computeShader.SetInt("resInt",resInt);
        
        rhoComputeBuffer = new ComputeBuffer(resInt*resInt, sizeof(float));
        fComputeBuffer =  new ComputeBuffer(resInt*resInt*4, sizeof(float));
        ftmpComputeBuffer =  new ComputeBuffer(resInt*resInt*4, sizeof(float));

        computeShader.SetBuffer(kernelInit, "rho", rhoComputeBuffer);
        computeShader.SetBuffer(kernelInit, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "ftmp", ftmpComputeBuffer);
        computeShader.SetBuffer(kernelStreaming, "ftmp", ftmpComputeBuffer);

        computeShader.SetBuffer(kernelCollision, "rho", rhoComputeBuffer);
        computeShader.SetBuffer(kernelCollision, "f", fComputeBuffer);

        computeShader.SetBuffer(kernelStreaming, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelStreaming, "rho", rhoComputeBuffer);

        computeShader.SetBuffer(kernelConstantBoundary, "f", fComputeBuffer);

        computeShader.SetBuffer(kernelGetRho, "f", fComputeBuffer);
        computeShader.SetBuffer(kernelGetRho, "rho", rhoComputeBuffer);

        threadGroups = (int) ((resInt*resInt + (threadGroupSize - 1)) / threadGroupSize);
        threadGroupsResInt = (int) ((resInt + (threadGroupSize - 1)) / threadGroupSize);
        
        computeShader.Dispatch(kernelInit, threadGroups, 1, 1);

        f = new float[resInt*resInt*4];
        ftmp = new float[resInt*resInt*4];
        rho = new float[resInt*resInt];

        rhoComputeBuffer.GetData(rho);
        UpdatePlot();
    }

    void LBMStep()
    {
        // Collision Process
        computeShader.Dispatch(kernelCollision, threadGroups, 1, 1);
        // Streaming Process
        computeShader.Dispatch(kernelStreaming, threadGroups, 1, 1);
        // Boundary Conditions
        computeShader.Dispatch(kernelConstantBoundary, threadGroupsResInt, 1, 1);
        // Get Rho
        computeShader.Dispatch(kernelGetRho, threadGroups, 1, 1);
    }
    
    void Update()
    {
        omega = 1f/(2f*alpha/(dt*csq)+0.5f);
        computeShader.SetFloat("omega",omega);
        computeShader.SetFloat("wallTemp1",wallTemp1);
        computeShader.SetFloat("wallTemp2",wallTemp2);
        computeShader.SetFloat("wallTemp3",wallTemp3);
        computeShader.SetFloat("wallTemp4",wallTemp4);
        
        for (int i = 0; i < loopCount; i++)
        {
            LBMStep();
        }
        UpdatePlot();
    }

    void UpdatePlot()
    {
        rhoComputeBuffer.GetData(rho);
        for (int i = 0; i < plotPixels.Length; i++)
        {
            plotPixels[i] = colorHeatMap.GetColorForValue(rho[i%resInt + (i/resInt)*resInt],1f);
        }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
    }
}
