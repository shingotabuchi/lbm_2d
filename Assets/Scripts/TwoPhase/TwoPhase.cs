// using System.Collections;
// using System;
// using System.Collections.Generic;
// using UnityEngine;
// using UnityEngine.UI;

// public class TwoPhase : MonoBehaviour
// {
//     public enum PlotMode
//     {
//         OrderParameter,
//         Density,
//         Speed,
//     }
//     public Image plotImage;
//     Texture2D plotTexture;
//     public int DIM_X = 50;
//     public int DIM_Y = 50;
//     int init,plotOrderParameter,collision,streaming;
//     ComputeBuffer uv,f,g,phi;
//     public ComputeShader compute;
//     RenderTexture renderTexture;
//     public int loopCount = 1;

//     public float maxRho,minRho;
//     public float maxSpeed,minSpeed;
//     public float maxPhi,minPhi;
//     public float tauf =  0.7f;
//     public float taug = 0.7f;
//     public float gamma = 10.0f, wid = 5.0f, phi0 = 1.0f;
//     public float sig = 0.0001f;
//     public float radius = 0.25f;

//     private void OnDestroy() {
//         uv.Dispose();
//         f.Dispose();
//         g.Dispose();
//         phi.Dispose();
//     }

//     private void Start() {
//         plotTexture = new Texture2D(DIM_X,DIM_Y);
//         plotTexture.filterMode = FilterMode.Point;
//         plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),UnityEngine.Vector2.zero);
//         ((RectTransform)plotImage.transform).sizeDelta = new Vector2((DIM_X*1080)/DIM_Y,1080);

//         renderTexture = new RenderTexture(DIM_X,DIM_Y,24);
//         renderTexture.enableRandomWrite = true;

//         uv = new ComputeBuffer(DIM_X*DIM_Y*2,sizeof(float));
//         f = new ComputeBuffer(DIM_X*DIM_Y*9*2,sizeof(float));
//         g = new ComputeBuffer(DIM_X*DIM_Y*9*2,sizeof(float));
//         phi = new ComputeBuffer(DIM_X*DIM_Y,sizeof(float));

//         compute.SetInt("DIM_X",DIM_X);
//         compute.SetInt("DIM_Y",DIM_Y);
//         compute.SetInt("DIMSqrd",DIM_X*DIM_Y);
//         compute.SetInt("DIMSqrd9",DIM_X*DIM_Y*9);
//         OnValidate();

//         init = compute.FindKernel("Init");
//         compute.SetBuffer(init,"uv",uv);
//         compute.SetBuffer(init,"f",f);
//         compute.SetBuffer(init,"g",g);
//         compute.SetBuffer(init,"phi",phi);

//         plotOrderParameter = compute.FindKernel("PlotOrderParameter");
//         compute.SetBuffer(plotOrderParameter,"phi",phi);
//         compute.SetTexture(plotOrderParameter,"renderTexture",renderTexture);

//         collision = compute.FindKernel("Collision");
//         compute.SetBuffer(collision,"uv",uv);
//         compute.SetBuffer(collision,"f",f);
//         compute.SetBuffer(collision,"g",g);
//         compute.SetBuffer(collision,"phi",phi);

//         streaming = compute.FindKernel("Streaming");
//         compute.SetBuffer(streaming,"f",f);
//         compute.SetBuffer(streaming,"g",g);

//         compute.Dispatch(init,(DIM_X+7)/8,(DIM_Y+7)/8,1);
//         compute.Dispatch(plotOrderParameter,(DIM_X+7)/8,(DIM_Y+7)/8,1);
//         RenderTexture.active = renderTexture;
//         plotTexture.ReadPixels(new Rect(0, 0, renderTexture.width, renderTexture.height), 0, 0);
//         plotTexture.Apply();
//     }

//     private void FixedUpdate() {
//         for (int kk = 0; kk < loopCount; kk++)
//         {
//             compute.Dispatch(collision,(DIM_X+7)/8,(DIM_Y+7)/8,1);
//             compute.Dispatch(streaming,(DIM_X+7)/8,(DIM_Y+7)/8,1);
//         }
//         compute.Dispatch(plotOrderParameter,(DIM_X+7)/8,(DIM_Y+7)/8,1);
//         RenderTexture.active = renderTexture;
//         plotTexture.ReadPixels(new Rect(0, 0, renderTexture.width, renderTexture.height), 0, 0);
//         plotTexture.Apply();
//     }

//     private void OnValidate() {
//         float beta  = 3.0f/4.0f*sig/wid*Mathf.Pow(phi0,4f);
//         float kap   = 3.0f/8.0f*sig*wid/Mathf.Pow(phi0,2f);

//         compute.SetFloat("maxRho",maxRho);
//         compute.SetFloat("minRho",minRho);
//         compute.SetFloat("maxSpeed",maxSpeed);
//         compute.SetFloat("minSpeed",minSpeed);
//         compute.SetFloat("maxPhi",maxPhi);
//         compute.SetFloat("minPhi",minPhi);

//         compute.SetFloat("beta",beta);
//         compute.SetFloat("kap",kap);
//         compute.SetFloat("tauf",tauf);
//         compute.SetFloat("taug",taug);
//         compute.SetFloat("radius",radius);
//         compute.SetFloat("phi0",phi0);
//         compute.SetFloat("gamma",gamma);
//         compute.SetFloat("wid",wid);
//     }
// }
