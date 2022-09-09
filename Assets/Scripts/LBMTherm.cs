using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.Numerics;
using System.Text;
using System.IO;
public class LBMTherm : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public HeatMapMode mode = HeatMapMode.Speed;
    public int DIM_X = 65;
    public int DIM_Y = 65;

    public float maxTemp,minTemp;
    float[,] e, ea, u, v;
    float[,,] g, g0, gtmp;
    float[] cx = new float[5]{0f,1f,0f,-1f, 0f};
    float[] cy = new float[5]{0f,0f,1f, 0f,-1f};
    float[] w = new float[5]{1f/3f,1f/6f,1f/6f,1f/6f,1f/6f};
    float h, tmp, norm, err, emax, emin,wav,chi,u0;
    public float pe = 20.0f, taug = 0.56f, q = 0.7f;
    public BoundaryTypeTherm boundaryType = BoundaryTypeTherm.Dirichlet;
    Complex beta;
    Complex[,] ec;
    
    public int loopCount = 1;

    // Start is called before the first frame update
    void Start()
    {
        plotTexture = new Texture2D(DIM_X,DIM_Y);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),UnityEngine.Vector2.zero);

        wav = 2f*Mathf.PI/(float)(DIM_X-1);
        chi = (taug - 0.5f)/3.0f;
        h = (float)(DIM_Y - 2) + 2f*q;
        u0 = pe*chi/h;
        beta = wav * Complex.Sqrt(new Complex(1.0f,u0/chi/wav));

        u = new float[DIM_X,DIM_Y];
        v = new float[DIM_X,DIM_Y];
        e = new float[DIM_X,DIM_Y];
        ea = new float[DIM_X,DIM_Y];
        ec = new Complex[DIM_X,DIM_Y];

        g = new float[5,DIM_X,DIM_Y];
        g0 = new float[5,DIM_X,DIM_Y];
        gtmp = new float[5,DIM_X,DIM_Y];

        maxTemp = 0f;
        minTemp = Mathf.Infinity;
        var sb = new StringBuilder("");
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u[i,j] = u0;
                v[i,j] = 0.0f;

                switch (boundaryType)
                {
                    case BoundaryTypeTherm.Dirichlet:
                    ec[i,j] = Complex.Exp(new Complex(0.0f,wav*(double)i))/(Complex.Exp(beta*h) - Complex.Exp(-beta*h))
                            *( (1.0f - Complex.Exp(-beta*h))*Complex.Exp( beta*((double)j - 1.0f + q))
                                - (1.0f - Complex.Exp( beta*h))*Complex.Exp(-beta*((double)j - 1.0f + q)));
                    break;
                    case BoundaryTypeTherm.Neumann:
                    ec[i,j] = Complex.Exp(new Complex(0.0f,wav*(double)i))/beta/h/(Complex.Exp(beta*h) - Complex.Exp(-beta*h))
                            *( (1.0f + Complex.Exp(-beta*h))*Complex.Exp( beta*((double)j - 1.0f + q))
                                + (1.0f + Complex.Exp( beta*h))*Complex.Exp(-beta*((double)j - 1.0f + q)));
                    break;
                }

                ea[i,j] = (float)ec[i,j].Real;
                e[i,j] = ea[i,j];
                print(e[i,j]);
                sb.Append(e[i,j].ToString("0.000000") + "\n");
                maxTemp = Mathf.Max(maxTemp,e[i,j]);
                minTemp = Mathf.Min(minTemp,e[i,j]);
                for (int k = 0; k < 5; k++)
                {
                    tmp = cx[k]*u[i,j] + cy[k]*v[i,j];      
                    if(k == 0) g0[k,i,j] = e[i,j]/3.0f;
                    else g0[k,i,j] = e[i,j]*(1.0f + 3.0f*tmp)/6.0f;
                    g[k,i,j] = g0[k,i,j];
                }
            }
        }
        // using (var writer = new StreamWriter("C:/lbm_2d/Assets/Scripts/data.txt", false))
        // {
        //     writer.Write(sb.ToString());
        // }
        UpdatePlot();
    }

    void LBMStep()
    {
        Collision();
        Streaming();
        HalfWayBounceBackBoundaries();
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
        for (int i = 0; i < plotPixels.Length; i++)
        {
            plotPixels[i] = colorHeatMap.GetColorForValue(e[i%DIM_X,i/DIM_X]-minTemp,maxTemp-minTemp);
        }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
    }

    void Collision()
    {
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                for (int k = 0; k < 5; k++)
                {
                    tmp = cx[k]*u[i,j] + cy[k]*v[i,j];      
                    if(k == 0) g0[k,i,j] = e[i,j]/3.0f;
                    else g0[k,i,j] = e[i,j]*(1.0f + 3.0f*tmp)/6.0f;
                    g[k,i,j] = g[k,i,j] - (g[k,i,j] - g0[k,i,j])/taug;
                }
            }   
        }
    }
    
    void Streaming()
    {
        gtmp = (float[,,])(g.Clone());
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            { 
                for(int k = 0; k < 5; k++)
                {
                    //periodic boundary condition
                    int im = (i + (int)cx[k] + DIM_X)%DIM_X; 
                    int jm = j + (int)cy[k];
                    if(jm!=DIM_Y&&jm!=-1)g[k,im,jm] = gtmp[k,i,j];
                } 
            }
        }
    }

    void HalfWayBounceBackBoundaries()
    {
        for (int i = 0; i < DIM_X; i++)
        {
            if(boundaryType == BoundaryTypeTherm.Dirichlet)
            {
                g[2,i,   1] = 2.0f*(q - 1.0f)*g[4,i, 0]
                      - Mathf.Pow((2.0f*q - 1.0f),2f)/(2.0f*q + 1.0f)*g[4,i, 1]
                      + 2.0f*(2.0f*q - 1.0f)   /(2.0f*q + 1.0f)*g[2,i, 2]
                      + (3.0f - 2.0f*q)/(2.0f*q + 1.0f)*Mathf.Cos(wav*(float)i)/3.0f;

                g[4,i,DIM_Y-2] = 2.0f*(    q - 1.0f)                 *g[2,i,DIM_Y-1  ]
                            - Mathf.Pow((2.0f*q - 1.0f),2)/(2.0f*q + 1.0f)*g[2,i,DIM_Y-2]
                            + 2.0f*(2.0f*q - 1.0f)   /(2.0f*q + 1.0f)*g[4,i,DIM_Y-3]
                            + (3.0f - 2.0f*q)/(2.0f*q + 1.0f)*Mathf.Cos(wav*(float)i)/3.0f;
            }
            else
            {
                g[2,i,   1] = g[4,i,0]
                            - (2.0f*q - 1.0f)/(2.0f*q + 1.0f)*g[4,i, 1]
                            + (2.0f*q - 1.0f)/(2.0f*q + 1.0f)*g[2,i, 2]
                            +           2.0f/(2.0f*q + 1.0f)*Mathf.Cos(wav*(float)i)/h*chi;
                g[4,i,DIM_Y-2] = g[2,i,DIM_Y-1]
                            - (2.0f*q - 1.0f)/(2.0f*q + 1.0f)*g[2,i,DIM_Y-2]
                            + (2.0f*q - 1.0f)/(2.0f*q + 1.0f)*g[4,i,DIM_Y-3]
                            +           2.0f/(2.0f*q + 1.0f)*Mathf.Cos(wav*(float)i)/h*chi;
            }
        }
    }

    void UpdateSpeedAndTemperature()
    {
        maxTemp = 0f;
        minTemp = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                e[i,j] = g[0,i,j]; 
                for(int k = 1; k <= 4; k++)
                {
                    e[i,j] = e[i,j] + g[k,i,j];
                    // u[i,j] =   u[i,j] + g[k,i,j]*cx[k];
                    // v[i,j] =   v[i,j] + g[k,i,j]*cy[k];
                } 
                // u[i,j] = u[i,j]/e[i,j];
                // v[i,j] = v[i,j]/e[i,j];

                maxTemp = Mathf.Max(maxTemp,e[i,j]);
                minTemp = Mathf.Min(minTemp,e[i,j]);
            } 
        }
    }
}
