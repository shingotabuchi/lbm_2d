using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LBMNatConvP1D : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public VectorField vectorField;
    public HeatMapMode mode = HeatMapMode.Speed;
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

    float[] cx = new float[9]{0, 1,    0,   -1,    0,     1,    -1,    -1,     1};
    float[] cy = new float[9]{0, 0,    1,    0,   -1,     1,     1,    -1,    -1};
    float[] wg = new float[5]{1f/3f,1f/6f,1f/6f,1f/6f,1f/6f};
    float[] wf = new float[9]{4f/9f,1f/9f,1f/9f,1f/9f,1f/9f,1f/36f,1f/36f,1f/36f,1f/36f};
    float[] rho, u, v, temp, fx, fy,speed,forceFromParticlesX,forceFromParticlesY;
    float[] f, f0, ftmp;
    float[] g, g0, gtmp;
    
    float umax, umin, tmp, u2, nu, chi, norm, taug, rbetag, h;
    public float pr = 0.71f;
    public float ra =   10000.0f;
    public float tauf = 0.8f;
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

    public int loopCount = 1;

    public bool noParticles = false;
    public bool particleConstantNoTemp = false;

    public Dropdown wallDropDown1;
    public Dropdown wallDropDown2;
    public Dropdown wallDropDown3;
    public Dropdown wallDropDown4;
    
    public Text wallText1;
    public Text wallText2;
    public Text wallText3;
    public Text wallText4;

    public Slider wallSlider1;
    public Slider wallSlider2;
    public Slider wallSlider3;
    public Slider wallSlider4;

    public Toggle applyChange;
    public Toggle normalize;

    public GameObject heatMap;

    Vector2 heatMapOriginalPos;
    float heatMapInitHeight;
    float heatMapHeight;

    public Transform canvas;

    public Dropdown heatMapModeDropDown;
    public Dropdown particleTempDropDown;
    public Dropdown initModeDropDown;
    public Slider particleTempSlider;
    public Slider particleRadiusSlider;

    public Toggle particleToggle;
    public Slider particleCountSlider;
    // Start is called before the first frame update
    void Start()
    {
        heatMapOriginalPos = heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().anchoredPosition;
        heatMapInitHeight = heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().sizeDelta.y;
        heatMapHeight = heatMapInitHeight;
        heatMap.transform.Find("MaxKnob").GetComponent<RectTransform>().anchoredPosition =
        new Vector2(
            heatMap.transform.Find("MaxKnob").GetComponent<RectTransform>().anchoredPosition.x,
            heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().anchoredPosition.y + 
            heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().sizeDelta.y/2f
        );

        heatMap.transform.Find("MinKnob").GetComponent<RectTransform>().anchoredPosition =
        new Vector2(
            heatMap.transform.Find("MinKnob").GetComponent<RectTransform>().anchoredPosition.x,
            heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().anchoredPosition.y - 
            heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().sizeDelta.y/2f
        );

        plotTexture = new Texture2D(DIM_X,DIM_Y);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),UnityEngine.Vector2.zero);

        if(DIM_X>DIM_Y)
        {
            plotImage.transform.GetComponent<RectTransform>().sizeDelta = new Vector2(
                plotImage.transform.GetComponent<RectTransform>().sizeDelta.x,
                (plotImage.transform.GetComponent<RectTransform>().sizeDelta.x*(float)DIM_Y)/(float)DIM_X
            );
        }
        if(DIM_X<DIM_Y)
        {
            plotImage.transform.GetComponent<RectTransform>().sizeDelta = new Vector2(
                (plotImage.transform.GetComponent<RectTransform>().sizeDelta.y*(float)DIM_X)/(float)DIM_Y,
                plotImage.transform.GetComponent<RectTransform>().sizeDelta.y
            );
        }

        vectorField.aspectRatio = (float)DIM_X/(float)DIM_Y;

        h = (float)(DIM_X-1 - 1);
        nu = (tauf - 0.5f)/3.0f;
        chi = nu/pr;
        taug = 3.0f*chi + 0.5f;

        rbetag = ra*nu*chi/h/h/h;
        speed = new float[DIM_X*DIM_Y];
        u = new float[DIM_X*DIM_Y];
        v = new float[DIM_X*DIM_Y];
        fx = new float[DIM_X*DIM_Y];
        forceFromParticlesX = new float[DIM_X*DIM_Y];
        forceFromParticlesY = new float[DIM_X*DIM_Y];
        fy = new float[DIM_X*DIM_Y];
        rho = new float[DIM_X*DIM_Y];
        temp = new float[DIM_X*DIM_Y];
        f = new float[9*DIM_X*DIM_Y];
        f0 = new float[9*DIM_X*DIM_Y];
        ftmp = new float[9*DIM_X*DIM_Y];
        g = new float[5*DIM_X*DIM_Y];
        g0 = new float[5*DIM_X*DIM_Y];
        gtmp = new float[5*DIM_X*DIM_Y];
        gravity[0] = 0f;
        gravity[1] = -rbetag/beta;
        
        float[] spawnPos = new float[particleCount*2];
        GameObject[] objs = new GameObject[particleCount];
        particleRadiusSlider.value = particleRadius;
        float scale = (2f*particleRadius*plotImage.transform.GetComponent<RectTransform>().sizeDelta.x * plotImage.transform.localScale.x)/(DIM_X*particlePrefab.GetComponent<RectTransform>().sizeDelta.x);
        plotSizeDelta = plotImage.transform.GetComponent<RectTransform>().sizeDelta;
        plotScale = plotImage.transform.localScale.x;
        particleCountSlider.value = particleCount;
        for (int i = 0; i < particleCount; i++)
        {
            spawnPos[i + 0*particleCount] = UnityEngine.Random.Range(particleRadius,DIM_X-particleRadius);
            spawnPos[i + 1*particleCount] = UnityEngine.Random.Range(particleRadius,DIM_Y-particleRadius);
            objs[i] = Instantiate(particlePrefab,particleParent);
            objs[i].GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
                                                            + new Vector2((spawnPos[i + 0*particleCount]*plotSizeDelta.x)/DIM_X,(spawnPos[i + 1*particleCount]*plotSizeDelta.y)/DIM_Y))*plotScale;
            objs[i].transform.localScale = new Vector3(scale,scale,1);
        }
        roundParticles = new RoundParticlesForCompute(DIM_X,particleCount,particleDensity,particleRadius,spawnPos,objs,Random.Range(0f,2f*Mathf.PI));
        
        maxTemp = 0f;
        minTemp = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u[i + j*DIM_X] = 0.0f; v[i + j*DIM_X] = 0.0f; 
                rho[i + j*DIM_X] = 1.0f;

                u2 = u[i + j*DIM_X]*u[i + j*DIM_X] + v[i + j*DIM_X]*v[i + j*DIM_X];   

                forceFromParticlesX[i + j*DIM_X] = 0f;
                forceFromParticlesY[i + j*DIM_X] = 0f;
                if(initModeDropDown.value == 0) temp[i + j*DIM_X] = 0f;
                else if(initModeDropDown.value == 1) temp[i + j*DIM_X] = (float)(DIM_X-1 - i)/(float)(DIM_X-1 - 1);
                else if(initModeDropDown.value == 2) temp[i + j*DIM_X] = (float)(DIM_Y-1 - j)/(float)(DIM_Y-1 - 1);
                
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i + j*DIM_X] + cy[k]*v[i + j*DIM_X];     
                    f0[k + (i + j*DIM_X)*9] = wf[k]*rho[i + j*DIM_X]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
                    f[k + (i + j*DIM_X)*9] = f0[k + (i + j*DIM_X)*9];
                    if(k<5)
                    {
                        g0[k + (i + j*DIM_X)*5] = wg[k]*temp[i + j*DIM_X]*(1.0f + 3.0f*tmp);
                        g[k + (i + j*DIM_X)*5] = g0[k + (i + j*DIM_X)*5];
                    }
                }
            }
        }
    }

    void LBMStep()
    {
        Collision();
        // Force();
        Streaming();
        Boundaries();
        UpdateSpeedAndTemperature();

        if(!noParticles)
        {
            ImmersedBoundary();
        }
    }

    // Update is called once per frame
    void Update()
    {
        normalizeHeatMap = normalize.isOn;
        if(heatMapModeDropDown.value == 0) mode = HeatMapMode.Temperature;
        else mode = HeatMapMode.Speed;

        if(applyChange.isOn)
        {
            if((int)particleCountSlider.value != particleCount)
            {
                int newparticleCount = (int)particleCountSlider.value;
                float[] spawnPos = new float[newparticleCount*2];
                float[] newTheta = new float[newparticleCount];
                GameObject[] objs = new GameObject[newparticleCount];
                float scale = (2f*particleRadius*plotImage.transform.GetComponent<RectTransform>().sizeDelta.x * plotImage.transform.localScale.x)/(DIM_X*particlePrefab.GetComponent<RectTransform>().sizeDelta.x);
                for (int i = 0; i < newparticleCount; i++)
                {
                    
                    if(i<particleCount)
                    {
                        newTheta[i] = roundParticles.theta[i];
                        spawnPos[i + 0*newparticleCount] = roundParticles.pos[i + 0*particleCount];
                        spawnPos[i + 1*newparticleCount] = roundParticles.pos[i + 1*particleCount];
                        objs[i] = roundParticles.objs[i];
                    }
                    else
                    {
                        newTheta[i] = Random.Range(0f,2f*Mathf.PI);
                        spawnPos[i + 0*newparticleCount] = UnityEngine.Random.Range(particleRadius,DIM_X-particleRadius);
                        spawnPos[i + 1*newparticleCount] = UnityEngine.Random.Range(particleRadius,DIM_Y-particleRadius);
                        objs[i] = Instantiate(particlePrefab,particleParent);
                        objs[i].GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
                                                                        + new Vector2((spawnPos[i + 0*newparticleCount]*plotSizeDelta.x)/DIM_X,(spawnPos[i + 1*newparticleCount]*plotSizeDelta.y)/DIM_Y))*plotScale;
                        objs[i].transform.localScale = new Vector3(scale,scale,1);
                    }
                }
                if(newparticleCount < particleCount)
                {
                    for (int i = newparticleCount; i < particleCount; i++)
                    {
                        Destroy(roundParticles.objs[i]);
                    }
                }
                particleCount = newparticleCount;
                roundParticles.UpdateParticleCount(particleCount,spawnPos,newTheta,objs);
            }
            if(particleRadiusSlider.value != particleRadius)
            {
                particleRadius = particleRadiusSlider.value;
                float scale = (2f*particleRadius*plotImage.transform.GetComponent<RectTransform>().sizeDelta.x * plotImage.transform.localScale.x)/(DIM_X*particlePrefab.GetComponent<RectTransform>().sizeDelta.x);
                for (int i = 0; i < particleCount; i++)
                {
                    roundParticles.objs[i].transform.localScale = new Vector3(scale,scale,1);
                }
                roundParticles.UpdateRadius(particleRadius);
            }
            if(particleToggle.isOn==noParticles)
            {
                noParticles = !particleToggle.isOn;
                particleParent.gameObject.SetActive(particleToggle.isOn);
            }
            if(particleTempDropDown.value == 0) particleConstantNoTemp = true;
            else particleConstantNoTemp = false;
            particleTemp = particleTempSlider.value;
            if(wallDropDown1.value == 0) wallboundaries[0] = BoundaryType.Adiabatic;
            else wallboundaries[0] = BoundaryType.Constant;
            if(wallDropDown2.value == 0) wallboundaries[1] = BoundaryType.Adiabatic;
            else wallboundaries[1] = BoundaryType.Constant;
            if(wallDropDown3.value == 0) wallboundaries[2] = BoundaryType.Adiabatic;
            else wallboundaries[2] = BoundaryType.Constant;
            if(wallDropDown4.value == 0) wallboundaries[3] = BoundaryType.Adiabatic;
            else wallboundaries[3] = BoundaryType.Constant;

            wallTemp1 = wallSlider1.value;
            wallTemp2 = wallSlider2.value;
            wallTemp3 = wallSlider3.value;
            wallTemp4 = wallSlider4.value;

            wallText1.text = "Wall temp : " + wallTemp1.ToString("0.000");
            wallText2.text = "Wall temp : " + wallTemp2.ToString("0.000");
            wallText3.text = "Wall temp : " + wallTemp3.ToString("0.000");
            wallText4.text = "Wall temp : " + wallTemp4.ToString("0.000");
        }

        for (int i = 0; i < loopCount; i++)
        {
            LBMStep();
        }
        UpdatePlot();
    }

    void UpdatePlot()
    {
        if(maxSpeed == 0f)maxSpeed = 1f;
        if(maxTemp == 0f)maxTemp = 1f;
        UpdateHeatMap();
        for (int i = 0; i < plotPixels.Length; i++)
        {
            if(mode == HeatMapMode.Speed) 
            plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X+(i/DIM_X)*DIM_X]-minSpeed,maxSpeed-minSpeed);
            else
            plotPixels[i] = colorHeatMap.GetColorForValue(temp[i%DIM_X+(i/DIM_X)*DIM_X]-minTemp,maxTemp-minTemp);
            // if(normalizeHeatMap)
            // {
                
            // }
            // else
            // {
            //     if(mode == HeatMapMode.Speed) 
            //     plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X+(i/DIM_X)*DIM_X],maxSpeed);
            //     else
            //     plotPixels[i] = colorHeatMap.GetColorForValue(temp[i%DIM_X+(i/DIM_X)*DIM_X],maxTemp);
            // }
        }
        if(!noParticles)
        {
            for (int i = 0; i < particleCount; i++)
            {
                // roundParticles.PlotParticlePerimeter(ref plotPixels, particleColor);
                // roundParticles.PlotParticleFill(ref plotPixels);
                roundParticles.objs[i].GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
                                                                                    + new Vector2((roundParticles.pos[i + 0*particleCount]*plotSizeDelta.x)/DIM_X,(roundParticles.pos[i + 1*particleCount]*plotSizeDelta.y)/DIM_Y))*plotScale;
                roundParticles.objs[i].transform.rotation = Quaternion.Euler(0,0,((roundParticles.theta[i]*180f)/Mathf.PI));
            }
        }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
        vectorField.UpdateVectors(u,v,DIM_Y,maxSpeed);
    }

    void Collision()
    {
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u2 = u[i + j*DIM_X]*u[i + j*DIM_X] + v[i + j*DIM_X]*v[i + j*DIM_X];   
                fx[i + j*DIM_X] = forceFromParticlesX[i + j*DIM_X]; fy[i + j*DIM_X] = rbetag*(temp[i + j*DIM_X] - 0.5f) + forceFromParticlesY[i + j*DIM_X];
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i + j*DIM_X] + cy[k]*v[i + j*DIM_X];     
                    f0[k + (i + j*DIM_X)*9] = wf[k]*rho[i + j*DIM_X]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);

                    f[k + (i + j*DIM_X)*9] = f[k + (i + j*DIM_X)*9] - (f[k + (i + j*DIM_X)*9] - f0[k + (i + j*DIM_X)*9])/tauf +  + 3f*wf[k]*(cx[k]*fx[i + j*DIM_X] + cy[k]*fy[i + j*DIM_X]);;
                    ftmp[k + (i + j*DIM_X)*9] = f[k + (i + j*DIM_X)*9];
                    if(k<5)
                    {
                        g0[k + (i + j*DIM_X)*5] = wg[k]*temp[i + j*DIM_X]*(1.0f + 3.0f*tmp);
                        g[k + (i + j*DIM_X)*5] = g[k + (i + j*DIM_X)*5] - (g[k + (i + j*DIM_X)*5] - g0[k + (i + j*DIM_X)*5])/taug;
                        gtmp[k + (i + j*DIM_X)*5] = g[k + (i + j*DIM_X)*5];
                    }
                }
                forceFromParticlesX[i + j*DIM_X] = 0f;
                forceFromParticlesY[i + j*DIM_X] = 0f;
            }   
        }
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

        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            { 
                for(int k = 0; k < 9; k++)
                {
                    int im = i + (int)cx[k]; 
                    int jm = j + (int)cy[k];
                    if((jm!=DIM_Y&&jm!=-1) && (im!=DIM_X&&im!=-1))
                    {
                        f[k + (im + jm*DIM_X)*9] = ftmp[k + (i + j*DIM_X)*9];
                        if(k<5)
                        g[k + (im + jm*DIM_X)*5] = gtmp[k + (i + j*DIM_X)*5];
                    }
                } 
            }
        }
    }

    void Boundaries()
    {
        for (int j = 0; j < DIM_Y; j++)
        // for (int j = 1; j < DIM_Y-1; j++)
        {
            f[1 + (0 + j*DIM_X)*9] = f[3 + (0 + j*DIM_X)*9];
            f[5 + (0 + j*DIM_X)*9] = f[7 + (0 + j*DIM_X)*9];
            f[8 + (0 + j*DIM_X)*9] = f[6 + (0 + j*DIM_X)*9];
            f[3 + (DIM_X-1 + j*DIM_X)*9] = f[1 + (DIM_X-1 + j*DIM_X)*9]; 
            f[7 + (DIM_X-1 + j*DIM_X)*9] = f[5 + (DIM_X-1 + j*DIM_X)*9]; 
            f[6 + (DIM_X-1 + j*DIM_X)*9] = f[8 + (DIM_X-1 + j*DIM_X)*9]; 
            
            switch (wallboundaries[0])
            {
                case BoundaryType.Constant: 
                g[0 + (0 + j*DIM_X)*5] = wg[0]*wallTemp1;
                g[1 + (0 + j*DIM_X)*5] = 2f*wg[1]*wallTemp1 - g[2 + (0 + j*DIM_X)*5];
                g[3 + (0 + j*DIM_X)*5] = 2f*wg[1]*wallTemp1 - g[4 + (0 + j*DIM_X)*5];
                break;
                case BoundaryType.Adiabatic:
                for (int k = 0; k < 5; k++)
                {
                    g[k + (0 + j*DIM_X)*5] = g[k+ (1 + j*DIM_X)*5];
                }
                break;
                case BoundaryType.Bounceback:
                g[1 + (0 + j*DIM_X)*5] = g[3 + (0 +j*DIM_X)*5 ];
                break;
            }
            switch (wallboundaries[1])
            {
                case BoundaryType.Constant: 
                g[0 + (DIM_X-1 + j*DIM_X)*5] = wg[0]*wallTemp2;
                g[2 + (DIM_X-1 + j*DIM_X)*5] = 2f*wg[1]*wallTemp2 - g[1 + (DIM_X-1 + j*DIM_X)*5];
                g[3 + (DIM_X-1 + j*DIM_X)*5] = 2f*wg[1]*wallTemp2 - g[4 + (DIM_X-1 + j*DIM_X)*5];
                break;
                case BoundaryType.Adiabatic:
                for (int k = 0; k < 5; k++)
                {
                    g[k + (DIM_X-1 + j*DIM_X)*5] = g[k+ (DIM_X-2 + j*DIM_X)*5];
                }
                break;
                case BoundaryType.Bounceback:
                g[3 + (DIM_X-1 + j*DIM_X)*5] =g[1+ (DIM_X-1 + j*DIM_X)*5];
                break;
            }
            // int jm = j - 1; int jp = j + 1;
            // f[1,    + (1 + j*DIM_X)*9] = f[3,0,j ];
            // f[5,    + (1 + j*DIM_X)*9] = f[7,0,jm];
            // f[8,    + (1 + j*DIM_X)*9] = f[6,0,jp];

            // f[3,DIM_X-1- + (1 + j*DIM_X)*9] = f[1+ (DIM_X-1 + j*DIM_X)*9];
            // f[7,DIM_X-1- + (1 + j*DIM_X)*9] = f[5,DIM_X-1,jp];
            // f[6,DIM_X-1- + (1 + j*DIM_X)*9] = f[8,DIM_X-1,jm];

            // g[1,    + (1 + j*DIM_X)*5] =-g[3, 0,j ] + 1.0f/3.0f;
            // g[3,DIM_X-1- + (1 + j*DIM_X)*5] =-g[1,DIM_X-1,j ];
        }
        for(int i = 0; i < DIM_X; i++){
        // for(int i = 1; i < DIM_X-1; i++){
            f[4+ (i + (DIM_Y-1)*DIM_X)*9] = f[2+ (i + (DIM_Y-1)*DIM_X)*9];
            f[7+ (i + (DIM_Y-1)*DIM_X)*9] = f[5+ (i + (DIM_Y-1)*DIM_X)*9];
            f[8+ (i + (DIM_Y-1)*DIM_X)*9] = f[6+ (i + (DIM_Y-1)*DIM_X)*9]; 
            f[2+ (i + 0*DIM_X)*9] = f[4+ (i + 0*DIM_X)*9]; 
            f[5+ (i + 0*DIM_X)*9] = f[7+ (i + 0*DIM_X)*9]; 
            f[6+ (i + 0*DIM_X)*9] = f[8+ (i + 0*DIM_X)*9]; 

            // int im = i - 1; int ip = i + 1;
            // f[2,i,   1] = f[4+ (i + 0*DIM_X)*9];
            // f[5,i,   1] = f[7,im, 0];
            // f[6,i,   1] = f[8,ip, 0];

            // f[4,i,DIM_Y-1-1] = f[2+ (i + (DIM_Y-1)*DIM_X)*9];
            // f[7,i,DIM_Y-1-1] = f[5,ip,DIM_Y-1];
            // f[8,i,DIM_Y-1-1] = f[6,im,DIM_Y-1];

            // g[2,i,   1] = g[4+ (i + 0*DIM_X)*5];
            // g[4,i,DIM_Y-1-1] = g[2+ (i + (DIM_Y-1)*DIM_X)*5];

            switch (wallboundaries[2])
            {
                case BoundaryType.Constant: 
                g[0+ (i + (DIM_Y-1)*DIM_X)*5] = wg[0]*wallTemp3;
                g[1+ (i + (DIM_Y-1)*DIM_X)*5] = 2f*wg[1]*wallTemp3 - g[2+ (i + (DIM_Y-1)*DIM_X)*5];
                g[4+ (i + (DIM_Y-1)*DIM_X)*5] = 2f*wg[1]*wallTemp3 - g[3+ (i + (DIM_Y-1)*DIM_X)*5];
                break;
                case BoundaryType.Adiabatic:
                for (int k = 0; k < 5; k++)
                {
                    g[k+ (i + (DIM_Y-1)*DIM_X)*5] = g[k+ (i + (DIM_Y-2)*DIM_X)*5];
                }
                break;
                case BoundaryType.Bounceback:
                g[4+ (i + (DIM_Y-1)*DIM_X)*5] = g[2+ (i + (DIM_Y-1)*DIM_X)*5];
                break;
            }
            switch (wallboundaries[3])
            {
                case BoundaryType.Constant: 
                g[0+ (i + 0*DIM_X)*5] = wg[0]*wallTemp4;
                g[1+ (i + 0*DIM_X)*5] = 2f*wg[1]*wallTemp4 - g[2+ (i + 0*DIM_X)*5];
                g[3+ (i + 0*DIM_X)*5] = 2f*wg[1]*wallTemp4 - g[4+ (i + 0*DIM_X)*5];
                break;
                case BoundaryType.Adiabatic:
                for (int k = 0; k < 5; k++)
                {
                    g[k+ (i + 0*DIM_X)*5] = g[k+ (i + 1*DIM_X)*5];
                }
                break;
                case BoundaryType.Bounceback:
                g[2+ (i + 0*DIM_X)*5] = g[4+ (i + 0*DIM_X)*5];
                break;
            }
        }
    }

    void UpdateSpeedAndTemperature()
    {
        maxSpeed = 0f;
        minSpeed = Mathf.Infinity;
        if(!noParticles&&!particleConstantNoTemp) maxTemp = particleTemp;
        else maxTemp = 0f;
        minTemp = Mathf.Infinity;
        maxRho = 0f;
        minRho = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u[i + j*DIM_X] = 0f; v[i + j*DIM_X] = 0f;
                rho[i + j*DIM_X] = f[0+(i+j*DIM_X)*9]; 
                temp[i + j*DIM_X] =  g[0+(i+j*DIM_X)*5];
                for(int k = 1; k <= 8; k++)
                {
                    rho[i + j*DIM_X] = rho[i + j*DIM_X] + f[k + (i + j*DIM_X)*9];
                    u[i + j*DIM_X] =   u[i + j*DIM_X] + f[k + (i + j*DIM_X)*9]*cx[k];
                    v[i + j*DIM_X] =   v[i + j*DIM_X] + f[k + (i + j*DIM_X)*9]*cy[k];
                    if(k<5){
                    temp[i + j*DIM_X] =   temp[i + j*DIM_X] + g[k + (i + j*DIM_X)*5];}
                    
                } 
                u[i + j*DIM_X] = u[i + j*DIM_X]/rho[i + j*DIM_X];
                v[i + j*DIM_X] = v[i + j*DIM_X]/rho[i + j*DIM_X];
                speed[i + j*DIM_X] = Mathf.Sqrt(u[i + j*DIM_X]*u[i + j*DIM_X] + v[i + j*DIM_X]*v[i + j*DIM_X]);
                
                maxTemp = Mathf.Max(maxTemp,temp[i + j*DIM_X]);
                minTemp = Mathf.Min(minTemp,temp[i + j*DIM_X]);
                maxSpeed = Mathf.Max(maxSpeed,speed[i + j*DIM_X]);
                minSpeed = Mathf.Min(minSpeed,speed[i + j*DIM_X]);
                maxRho = Mathf.Max(maxRho,rho[i + j*DIM_X]);
                minRho = Mathf.Min(minRho,rho[i + j*DIM_X]);
            } 
        }
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
                float dist = roundParticles.ParticleDistance(n,k);
                // print(dist);
                for (int i = 0; i < 2; i++)
                {
                    if(dist < 2.0f*roundParticles.radius[n] + zeta)
                    {
                        roundParticles.forceFromCollisions[n + i*particleCount] += (roundParticles.pos[n + i*particleCount] - roundParticles.pos[k + i*particleCount])*(2.0f*roundParticles.radius[n] - dist + zeta)*(2.0f*roundParticles.radius[n] - dist + zeta)/epsw;
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
                if(!particleConstantNoTemp)
                {
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
                            forceFromParticlesX[i + j*DIM_X] += roundParticles.forceOnPerimeter[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] * tmp3 * 2.0f*Mathf.PI*roundParticles.radius[n]/(float)roundParticles.perimeterPointCount[n];
                            forceFromParticlesY[i + j*DIM_X] += roundParticles.forceOnPerimeter[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] * tmp3 * 2.0f*Mathf.PI*roundParticles.radius[n]/(float)roundParticles.perimeterPointCount[n];
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
            roundParticles.UpdatePosVel(n,gravity);
            roundParticles.UpdateOmegaTheta(n);
            roundParticles.UpdatePerimeter(n);
        }
    }

    void UpdateHeatMap()
    {
        if(mode == HeatMapMode.Temperature)
        {
            if(!normalizeHeatMap)
            {
                if(heatMap.transform.Find("MaxKnob").GetComponent<Knob>().pressed)
                {
                    heatMap.transform.Find("MaxKnob").GetComponent<RectTransform>().anchoredPosition = new Vector2(
                        heatMap.transform.Find("MaxKnob").GetComponent<RectTransform>().anchoredPosition.x,
                        MouseScreenPosition().y
                    );
                }
                if(heatMap.transform.Find("MinKnob").GetComponent<Knob>().pressed)
                {
                    heatMap.transform.Find("MinKnob").GetComponent<RectTransform>().anchoredPosition = new Vector2(
                        heatMap.transform.Find("MinKnob").GetComponent<RectTransform>().anchoredPosition.x,
                        MouseScreenPosition().y
                    );
                }
                float knobMax = (heatMap.transform.Find("MaxKnob").GetComponent<RectTransform>().anchoredPosition.y - heatMapOriginalPos.y + heatMapInitHeight/2f)/heatMapInitHeight;
                float knobMin = (heatMap.transform.Find("MinKnob").GetComponent<RectTransform>().anchoredPosition.y - heatMapOriginalPos.y + heatMapInitHeight/2f)/heatMapInitHeight;
                if(knobMax>maxTemp) maxTemp = knobMax;
                if(knobMin<minTemp) minTemp = knobMin;
            }
            float max = maxTemp;
            if(max>1f) max = 1f;
            if(minTemp<0f) minTemp = 0f;
            heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().sizeDelta = new Vector2(
                heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().sizeDelta.x,
                heatMapInitHeight*(max-minTemp)
            );
            heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().anchoredPosition =
            heatMapOriginalPos + 
            new Vector2(
                0,
                heatMapInitHeight*((max+minTemp)/2f - 0.5f)
            );

            heatMap.transform.Find("MaxKnob").GetComponent<RectTransform>().anchoredPosition =
            new Vector2(
                heatMap.transform.Find("MaxKnob").GetComponent<RectTransform>().anchoredPosition.x,
                heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().anchoredPosition.y + 
                heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().sizeDelta.y/2f
            );

            heatMap.transform.Find("MinKnob").GetComponent<RectTransform>().anchoredPosition =
            new Vector2(
                heatMap.transform.Find("MinKnob").GetComponent<RectTransform>().anchoredPosition.x,
                heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().anchoredPosition.y - 
                heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().sizeDelta.y/2f
            );

            heatMap.transform.Find("White").GetComponent<RectTransform>().sizeDelta = new Vector2(
                heatMap.transform.Find("White").GetComponent<RectTransform>().sizeDelta.x,
                heatMapInitHeight*(1f-(max+minTemp)/2f)
            );
            heatMap.transform.Find("White").GetComponent<RectTransform>().anchoredPosition =
            heatMapOriginalPos + 
            new Vector2(
                0,
                heatMapInitHeight*((max+minTemp)/4f)
            );
        }
        else
        {
            if(!normalizeHeatMap)
            {
                if(heatMap.transform.Find("MaxKnob").GetComponent<Knob>().pressed)
                {
                    heatMap.transform.Find("MaxKnob").GetComponent<RectTransform>().anchoredPosition = new Vector2(
                        heatMap.transform.Find("MaxKnob").GetComponent<RectTransform>().anchoredPosition.x,
                        MouseScreenPosition().y
                    );
                }
                if(heatMap.transform.Find("MinKnob").GetComponent<Knob>().pressed)
                {
                    heatMap.transform.Find("MinKnob").GetComponent<RectTransform>().anchoredPosition = new Vector2(
                        heatMap.transform.Find("MinKnob").GetComponent<RectTransform>().anchoredPosition.x,
                        MouseScreenPosition().y
                    );
                }
                float knobMax = (heatMap.transform.Find("MaxKnob").GetComponent<RectTransform>().anchoredPosition.y - heatMapOriginalPos.y + heatMapInitHeight/2f)/heatMapInitHeight;
                float knobMin = (heatMap.transform.Find("MinKnob").GetComponent<RectTransform>().anchoredPosition.y - heatMapOriginalPos.y + heatMapInitHeight/2f)/heatMapInitHeight;
                if(knobMax>maxSpeed) maxSpeed = knobMax;
                if(knobMin<minSpeed) minSpeed = knobMin;
            }
            float max = maxSpeed;
            if(max>1f) max = 1f;
            if(minSpeed<0f) minSpeed = 0f;
            heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().sizeDelta = new Vector2(
                heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().sizeDelta.x,
                heatMapInitHeight*(max-minSpeed)
            );
            heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().anchoredPosition =
            heatMapOriginalPos + 
            new Vector2(
                0,
                heatMapInitHeight*((max+minSpeed)/2f - 0.5f)
            );

            heatMap.transform.Find("MaxKnob").GetComponent<RectTransform>().anchoredPosition =
            new Vector2(
                heatMap.transform.Find("MaxKnob").GetComponent<RectTransform>().anchoredPosition.x,
                heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().anchoredPosition.y + 
                heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().sizeDelta.y/2f
            );

            heatMap.transform.Find("MinKnob").GetComponent<RectTransform>().anchoredPosition =
            new Vector2(
                heatMap.transform.Find("MinKnob").GetComponent<RectTransform>().anchoredPosition.x,
                heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().anchoredPosition.y - 
                heatMap.transform.Find("HeatMap").GetComponent<RectTransform>().sizeDelta.y/2f
            );

            heatMap.transform.Find("White").GetComponent<RectTransform>().sizeDelta = new Vector2(
                heatMap.transform.Find("White").GetComponent<RectTransform>().sizeDelta.x,
                heatMapInitHeight*(1f-(max+minSpeed)/2f)
            );
            heatMap.transform.Find("White").GetComponent<RectTransform>().anchoredPosition =
            heatMapOriginalPos + 
            new Vector2(
                0,
                heatMapInitHeight*((max+minSpeed)/4f)
            );
        }
    }

    public Vector2 MouseScreenPosition()
    {
        return new Vector2(
            (Input.mousePosition.x - canvas.GetComponent<RectTransform>().sizeDelta.x/2f),
            (Input.mousePosition.y - canvas.GetComponent<RectTransform>().sizeDelta.y/2f)
        );  
    }

    public void Reset()
    {
        float[] spawnPos = new float[particleCount*2];
        GameObject[] objs = new GameObject[particleCount];
        particleRadiusSlider.value = particleRadius;
        float scale = (2f*particleRadius*plotImage.transform.GetComponent<RectTransform>().sizeDelta.x * plotImage.transform.localScale.x)/(DIM_X*particlePrefab.GetComponent<RectTransform>().sizeDelta.x);
        plotSizeDelta = plotImage.transform.GetComponent<RectTransform>().sizeDelta;
        plotScale = plotImage.transform.localScale.x;
        particleCountSlider.value = particleCount;
        for (int i = 0; i < particleCount; i++)
        {
            Destroy(roundParticles.objs[i]);
            spawnPos[i + 0*particleCount] = UnityEngine.Random.Range(particleRadius,DIM_X-particleRadius);
            spawnPos[i + 1*particleCount] = UnityEngine.Random.Range(particleRadius,DIM_Y-particleRadius);
            objs[i] = Instantiate(particlePrefab,particleParent);
            objs[i].GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
                                                            + new Vector2((spawnPos[i + 0*particleCount]*plotSizeDelta.x)/DIM_X,(spawnPos[i + 1*particleCount]*plotSizeDelta.y)/DIM_Y))*plotScale;
            objs[i].transform.localScale = new Vector3(scale,scale,1);
        }
        roundParticles = new RoundParticlesForCompute(DIM_X,particleCount,particleDensity,particleRadius,spawnPos,objs,Random.Range(0f,2f*Mathf.PI));
        
        maxTemp = 0f;
        minTemp = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u[i + j*DIM_X] = 0.0f; v[i + j*DIM_X] = 0.0f; 
                rho[i + j*DIM_X] = 1.0f;

                u2 = u[i + j*DIM_X]*u[i + j*DIM_X] + v[i + j*DIM_X]*v[i + j*DIM_X];   

                forceFromParticlesX[i + j*DIM_X] = 0f;
                forceFromParticlesY[i + j*DIM_X] = 0f;
                if(initModeDropDown.value == 0) temp[i + j*DIM_X] = 0f;
                else if(initModeDropDown.value == 1) temp[i + j*DIM_X] = (float)(DIM_X-1 - i)/(float)(DIM_X-1 - 1);
                else if(initModeDropDown.value == 2) temp[i + j*DIM_X] = (float)(DIM_Y-1 - j)/(float)(DIM_Y-1 - 1);
                
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i + j*DIM_X] + cy[k]*v[i + j*DIM_X];     
                    f0[k + (i + j*DIM_X)*9] = wf[k]*rho[i + j*DIM_X]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
                    f[k + (i + j*DIM_X)*9] = f0[k + (i + j*DIM_X)*9];
                    if(k<5)
                    {
                        g0[k + (i + j*DIM_X)*5] = wg[k]*temp[i + j*DIM_X]*(1.0f + 3.0f*tmp);
                        g[k + (i + j*DIM_X)*5] = g0[k + (i + j*DIM_X)*5];
                    }
                }
            }
        }
    }
}