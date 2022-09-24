using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LBMCage : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public VectorField vectorField;
    public int DIM = 64;
    public HeatMapMode mode = HeatMapMode.Speed;

    public float tau = 0.56f;
    float tmp,u2,h,nu,norm,nb,maxrho,prevMaxSpeed;
    public float maxSpeed;
    float[] cx = new float[9]{0,        1,    0,   -1,    0,     1,    -1,    -1,     1};
    float[] cy = new float[9]{0,        0,    1,    0,   -1,     1,     1,    -1,    -1};
    float[] w = new  float[9]{4f/9f,1f/9f,1f/9f,1f/9f,1f/9f,1f/36f,1f/36f,1f/36f,1f/36f};
    float[,] u,v,rho,speed,gx,gy;
    float[,,] f,f0,ftmp;

    public int loopCount = 1;

    public Transform canvas;
    public float force;
    Vector3 prevMousePos;
    // Start is called before the first frame update
    void Start()
    {
        plotTexture = new Texture2D(DIM,DIM);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM,DIM),Vector2.zero);

        gx = new float[DIM,DIM];
        gy = new float[DIM,DIM];
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
        prevMousePos = MousePosition();
    }

    void LBMStep()
    {
        Collision();
        Streaming();
        BounceBackBoundaries();
        UpdateSpeedAndDensity();
    }

    // Update is called once per frame
    void Update()
    {
        CalculateForce();
        for (int i = 0; i < loopCount; i++)
        {
            LBMStep();
        }
        UpdatePlot();
    }

    void CalculateForce()
    {
        for(int i = 0; i < DIM; i++)
        { 
            for(int j = 0; j < DIM; j++)
            {
                gx[i,j] = 0f;
                gy[i,j] = 0f;
            }   
        }
        Vector3 currentMousePos = MousePosition();
        if(Input.GetMouseButton(0))
        {   
            int[] delta = new int[2] { ((int)currentMousePos.x - (int)prevMousePos.x), ((int)currentMousePos.y - (int)prevMousePos.y) };
            float xdiff,ydiff;
            xdiff = currentMousePos.x - prevMousePos.x;
            ydiff = currentMousePos.y - prevMousePos.y;
            if (delta[0] == 0)
            {
                for (int i = 0; i < (int)Mathf.Abs(delta[1]); i++)
                {
                    if (delta[1] > 0)
                    {
                        int im = (int)prevMousePos.x;
                        int jm = i + (int)prevMousePos.y;
                        if(im<0 || im>=DIM || jm<0 || jm>=DIM)continue;
                        gx[im,jm] = force*xdiff/Mathf.Sqrt(xdiff*xdiff + ydiff*ydiff);
                        gy[im,jm] = force*ydiff/Mathf.Sqrt(xdiff*xdiff + ydiff*ydiff);
                    }
                    else
                    {
                        int im = (int)prevMousePos.x;
                        int jm = -i + (int)prevMousePos.y;
                        if(im<0 || im>=DIM || jm<0 || jm>=DIM)continue;
                        gx[im,jm] = force*xdiff/Mathf.Sqrt(xdiff*xdiff + ydiff*ydiff);
                        gy[im,jm] = force*ydiff/Mathf.Sqrt(xdiff*xdiff + ydiff*ydiff);
                    }
                }
            }
            else if (delta[1] == 0)
            {
                for (int i = 0; i < (int)Mathf.Abs(delta[0]); i++)
                {
                    if (delta[0] > 0)
                    {
                        int im = i + (int)prevMousePos.x;
                        int jm = (int)prevMousePos.y;
                        if(im<0 || im>=DIM || jm<0 || jm>=DIM)continue;
                        gx[im,jm] = force*xdiff/Mathf.Sqrt(xdiff*xdiff + ydiff*ydiff);
                        gy[im,jm] = force*ydiff/Mathf.Sqrt(xdiff*xdiff + ydiff*ydiff);
                    }
                    else
                    {
                        int im = -i + (int)prevMousePos.x;
                        int jm = (int)prevMousePos.y;
                        if(im<0 || im>=DIM || jm<0 || jm>=DIM)continue;
                        gx[im,jm] = force*xdiff/Mathf.Sqrt(xdiff*xdiff + ydiff*ydiff);
                        gy[im,jm] = force*ydiff/Mathf.Sqrt(xdiff*xdiff + ydiff*ydiff);
                    }
                }
            }
            else
            {
                int x, y, dx, dy;
                x = 0;
                y = 0;
                if (delta[0] > 0) dx = 1;
                else dx = -1;
                if (delta[1] > 0) dy = 1;
                else dy = -1;
                while (x + (int)prevMousePos.x != (int)currentMousePos.x || y + (int)prevMousePos.y != (int)currentMousePos.y)
                {
                    int gaiseki_dx = (int)Mathf.Abs((x + dx) * delta[1] - y * delta[0]);
                    int gaiseki_dy = (int)Mathf.Abs(x * delta[1] - (y + dy) * delta[0]);

                    if (gaiseki_dx > gaiseki_dy)
                    {
                        y += dy;
                    }
                    else x += dx;

                    int im = x + (int)prevMousePos.x;
                    int jm = y + (int)prevMousePos.y;
                    if(im<0 || im>=DIM || jm<0 || jm>=DIM)continue;
                    gx[im,jm] = force*xdiff/Mathf.Sqrt(xdiff*xdiff + ydiff*ydiff);
                    gy[im,jm] = force*ydiff/Mathf.Sqrt(xdiff*xdiff + ydiff*ydiff);
                }
            }
        }
        prevMousePos = currentMousePos;
    }

    void UpdatePlot()
    {
        if(maxSpeed == 0f) maxSpeed = 1f;
        // else Debug.Log(maxSpeed);
        for (int i = 0; i < plotPixels.Length; i++)
        {
            // float c = speed[i%DIM,i/DIM]/maxSpeed;
            // plotPixels[i] = new Color(c,c,c,1);
            // if(speed[i%DIM,i/DIM]!=0f)Debug.Log(speed[i%DIM,i/DIM]);
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
                    f[k,i,j] += 3f*w[k]*rho[i,j]*(cx[k]*gx[i,j] + cy[k]*gy[i,j]);
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
                    // //periodic boundary condition
                    // int im = (i + (int)cx[k] + DIM)%DIM; 
                    // int jm = j + (int)cy[k];
                    // if(jm != DIM && jm !=-1) f[k,im,jm] = ftmp[k,i,j];
                    int im = i + (int)cx[k]; 
                    int jm = j + (int)cy[k];
                    if((jm!=DIM&&jm!=-1) && (im!=DIM&&im!=-1))
                    {
                        f[k,im,jm] = ftmp[k,i,j];
                    }
                } 
            }
        }
    }

    void BounceBackBoundaries()
    {
        // ftmp = (float[,,])(f.Clone());
        for (int i = 0; i < DIM; i++)
        {
            f[4,i,DIM-1] = f[2,i,DIM-1];
            f[7,i,DIM-1] = f[5,i,DIM-1];
            f[8,i,DIM-1] = f[6,i,DIM-1]; 
            f[2,i,0] = f[4,i,0]; 
            f[5,i,0] = f[7,i,0]; 
            f[6,i,0] = f[8,i,0]; 
        }
        for (int j = 0; j < DIM; j++)
        {
            f[3,DIM-1,j] = f[1,DIM-1,j];
            f[6,DIM-1,j] = f[8,DIM-1,j];
            f[7,DIM-1,j] = f[5,DIM-1,j];
            f[1,0,j] = f[3,0,j]; 
            f[5,0,j] = f[7,0,j]; 
            f[8,0,j] = f[6,0,j]; 
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

                speed[i,j] = Mathf.Sqrt(u[i,j]*u[i,j] + v[i,j]*v[i,j]);
                maxSpeed = Mathf.Max(maxSpeed,speed[i,j]);
                maxrho = Mathf.Max(maxrho,rho[i,j]);
            } 
        }
    }

    Vector3 MousePosition()
    {
        Vector2 plotSizeDelta = plotImage.transform.GetComponent<RectTransform>().sizeDelta;
        Vector3 plotScale = plotImage.transform.localScale;
        return new Vector3(
            ((((Input.mousePosition.x - canvas.GetComponent<RectTransform>().sizeDelta.x/2f) * canvas.localScale.x) + plotScale.x * plotSizeDelta.x/2f)*DIM)/(plotScale.x*plotSizeDelta.x),
            ((((Input.mousePosition.y - canvas.GetComponent<RectTransform>().sizeDelta.y/2f) * canvas.localScale.y) + plotScale.y * plotSizeDelta.y/2f)*DIM)/(plotScale.y*plotSizeDelta.y),
            0
        );  
    }
}
