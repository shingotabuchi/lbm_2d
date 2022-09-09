using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class FDMMain : MonoBehaviour
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

    float[,] tempCopy,temp;
    float dx = 1f;
    float dy = 1f;
    float dd;
    public int inversedt = 5;
    float dt;
    
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    void Start()
    {
        dt = 1f/(float)inversedt;
        dd = 1f/(dx*dx) + 1f/(dy*dy);
        plotTexture = new Texture2D(resInt,resInt);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,resInt,resInt),Vector2.zero);
        
        temp = new float[resInt,resInt];
        tempCopy = new float[resInt,resInt];
        for (int i = 0; i < resInt; i++)
        {
            for (int j = 0; j < resInt; j++)
            {
                if(i==0) temp[i,j] = wallTemp1;
                else temp[i,j] = 0f;
            }
        }
        UpdatePlot();
    }
    
    void Update()
    {
        for (int k = 0; k < (int)inversedt; k++)
        {
            tempCopy = (float[,])(temp.Clone());
            for (int i = 1; i < resInt-1; i++)
            {
                for (int j = 1; j < resInt-1; j++)
                {
                    float termx = (tempCopy[i+1,j] + tempCopy[i-1,j])/(dx*dx);
                    float termy = (tempCopy[i,j+1] + tempCopy[i,j-1])/(dy*dy);
                    temp[i,j] = tempCopy[i,j] + dt*alpha*(termx + termy - 2f*tempCopy[i,j]*dd);
                }
            }
            // Boundary Conditions
            for (int i = 0; i < resInt; i++)
            {
                temp[i,0] = temp[i,1];
            }
        }
        UpdatePlot();
    }

    void UpdatePlot()
    {
        for (int i = 0; i < plotPixels.Length; i++)
        {
            plotPixels[i] = colorHeatMap.GetColorForValue(temp[i%resInt,i/resInt],1f);
        }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
    }
}
