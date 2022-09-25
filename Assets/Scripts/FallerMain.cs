using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.SceneManagement;
public class FallerMain : MonoBehaviour
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
    Vector2 plotAnchPos; 
    Vector2 plotSizeDelta; 
    Vector3 plotScale; 
    [Range(0.0f, 1.0f)]
    public float particleTemp = 1f;

    // public int loopCount = 1;

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

    public bool isFaller = false;

    public int lineRes = 50;
    public int modes = 4;
    public float[] modeCoeffs = new float[5];
    public float[] modeSinCoeffs = new float[5];
    public float[] originalModeCoeffs = new float[5];
    public float[] originalModeSinCoeffs = new float[5];
    public float[] destModeCoeffs = new float[5];
    public float[] destModeSinCoeffs = new float[5];
    public float area = 1000f;
    public float realArea;
    public Slider[] CosSliders;
    public Slider[] SinSliders;
    public Slider[] CosShuukiSliders;
    public Slider[] SinShuukiSliders;
    public Color perimColor;
    public Color traceColor;
    public int time = 0;
    public float gx = 0f;
    public float gy = 0f;

    PolygonParticle polygonParticle;

    List<Vector3> dotPositions = new List<Vector3>();

    // public LineRenderer line;

    public float pTheta;

    public GameObject meshObj;
    public float segmentDivisor = 2f;
    Mesh mesh;
    public float coeffChangeSpeed = 0.001f;
    public Toggle shindoAll;
    public Slider sizeSlider;
    public InputField inputFieldH;
    public InputField inputFieldW;
    public Slider loopCountSlider;
    public Slider gxSlider;
    public Slider gySlider;
    public bool vicsekLike = false;

    public Toggle align;
    public Slider viewRangeSlider;
    public Slider noiseAngleSlider;
    public Slider particleFovSlider;
    public Slider particlePropelSlider;
    float particleViewRange,particleFov,noiseAngleDeg,particlePropulsion;
    
    public LineRenderer traceLine;
    List<Vector3> tracePoints = new List<Vector3>();
    public float traceThres = 0.5f;
    Vector2 originalPlotSizeDelta;
    public Toggle showTrace;
    // Start is called before the first frame update
    void Start()
    {
        if(isFaller)
        {
            mesh = new Mesh();
            meshObj.GetComponent<MeshFilter>().mesh = mesh;
        }
        if(vicsekLike)
        {
            particleViewRange = viewRangeSlider.value;
            particleFov = particleFovSlider.value;
            noiseAngleDeg = noiseAngleSlider.value;
            particlePropulsion = particlePropelSlider.value;
        }
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

        DIM_X = int.Parse(inputFieldW.text);
        DIM_Y = int.Parse(inputFieldH.text);
        plotTexture = new Texture2D(DIM_X,DIM_Y);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),UnityEngine.Vector2.zero);
        originalPlotSizeDelta = plotImage.transform.GetComponent<RectTransform>().sizeDelta;
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
        plotAnchPos = plotImage.transform.GetComponent<RectTransform>().anchoredPosition;
        plotSizeDelta = plotImage.transform.GetComponent<RectTransform>().sizeDelta;
        plotScale = plotImage.transform.localScale;
        if(!isFaller)
        {
            particleCount = (int)particleCountSlider.value;
            float[] spawnPos = new float[particleCount*2];
            GameObject[] objs = new GameObject[particleCount];
            particleRadiusSlider.value = particleRadius;
            float scale = (2f*particleRadius*plotImage.transform.GetComponent<RectTransform>().sizeDelta.x * plotImage.transform.localScale.x)/(DIM_X*particlePrefab.GetComponent<RectTransform>().sizeDelta.x);
            for (int i = 0; i < particleCount; i++)
            {
                spawnPos[i + 0*particleCount] = UnityEngine.Random.Range(particleRadius,DIM_X-particleRadius);
                spawnPos[i + 1*particleCount] = UnityEngine.Random.Range(particleRadius,DIM_Y-particleRadius);
                objs[i] = Instantiate(particlePrefab,particleParent);
                objs[i].GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
                                                                + new Vector2((spawnPos[i + 0*particleCount]*plotSizeDelta.x)/DIM_X,(spawnPos[i + 1*particleCount]*plotSizeDelta.y)/DIM_Y))*plotScale.x;
                objs[i].transform.localScale = new Vector3(scale,scale,1);
            }
            roundParticles = new RoundParticlesForCompute(DIM_X,particleCount,particleDensity,particleRadius,spawnPos,objs,Random.Range(0f,2f*Mathf.PI));
        }
        else
        {
            for (int i = 0; i < 4; i++)
            {
                if(i>=1) modeCoeffs[i] = CosSliders[i-1].value;
                destModeCoeffs[i] = modeCoeffs[i];
                originalModeCoeffs[i] = modeCoeffs[i];
                if(i>=1) modeSinCoeffs[i] = SinSliders[i-1].value;
                destModeSinCoeffs[i] = modeSinCoeffs[i];
                originalModeSinCoeffs[i]  = modeSinCoeffs[i];
            }
            UpdatePolygonPerimeter(true);
        }
        
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

        if(!noParticles&&!isFaller)
        {
            ImmersedBoundary();
        }
        if(isFaller)
        {
            FallerImmersedBoundary();

            for (int i = 0; i < 4; i++)
            {
                float diff = originalModeCoeffs[i]-destModeCoeffs[i];
                if(Mathf.Abs(diff) > coeffChangeSpeed)
                {
                    originalModeCoeffs[i] += (-diff/Mathf.Abs(diff))*coeffChangeSpeed;
                }
                else originalModeCoeffs[i] = destModeCoeffs[i];
                diff = originalModeSinCoeffs[i]-destModeSinCoeffs[i];
                if(Mathf.Abs(diff) > coeffChangeSpeed)
                {
                    originalModeSinCoeffs[i] += (-diff/Mathf.Abs(diff))*coeffChangeSpeed;
                }
                else originalModeSinCoeffs[i] = destModeSinCoeffs[i];
            }
            int shuukiFactor = 1;
            for (int i = 1; i < 4; i++)
            {
                shuukiFactor *= (int)CosShuukiSliders[i-1].value * (int)SinShuukiSliders[i-1].value;
                if(CosShuukiSliders[i-1].transform.Find("Toggle").GetComponent<Toggle>().isOn)
                {
                    modeCoeffs[i] = originalModeCoeffs[i] * Mathf.Cos(2f*Mathf.PI*((float)time/CosShuukiSliders[i-1].value));
                }
                else
                {   
                    modeCoeffs[i] = originalModeCoeffs[i];
                }
                if(SinShuukiSliders[i-1].transform.Find("Toggle").GetComponent<Toggle>().isOn)
                {
                    modeSinCoeffs[i] = originalModeSinCoeffs[i] * Mathf.Cos(2f*Mathf.PI*((float)time/SinShuukiSliders[i-1].value));
                }
                else
                {
                    modeSinCoeffs[i] = originalModeSinCoeffs[i];
                }
            }
            time = (time + 1)%shuukiFactor;
        }
    }

    // Update is called once per frame
    void Update()
    {
        normalizeHeatMap = normalize.isOn;
        if(heatMapModeDropDown.value == 0) mode = HeatMapMode.Temperature;
        else mode = HeatMapMode.Speed;

        if(isFaller)
        {
            traceLine.transform.gameObject.SetActive(showTrace.isOn);
            gxSlider.transform.Find("Text").GetComponent<Text>().text = "Gx=" + gxSlider.value.ToString("0.000");
            gySlider.transform.Find("Text").GetComponent<Text>().text = "Gy=" + gySlider.value.ToString("0.000");
            for (int i = 0; i < 3; i++)
            {
                CosSliders[i].transform.Find("Text").GetComponent<Text>().text = "a" + (i+1).ToString() + "= " + CosSliders[i].value.ToString("0.000");
                SinSliders[i].transform.Find("Text").GetComponent<Text>().text = "b" + (i+1).ToString() + "= " + SinSliders[i].value.ToString("0.000");
                CosShuukiSliders[i].transform.Find("Text").GetComponent<Text>().text = "period(step)= " + CosShuukiSliders[i].value.ToString("0");
                SinShuukiSliders[i].transform.Find("Text").GetComponent<Text>().text = "period(step)= " + SinShuukiSliders[i].value.ToString("0");
            }
        }
        else
        {
            particleCountSlider.transform.Find("Text").GetComponent<Text>().text = "Particle count = " + particleCountSlider.value.ToString("0");
            particleTempSlider.transform.Find("Text").GetComponent<Text>().text = "Particle temp = " + particleTempSlider.value.ToString("0.000");
            particleRadiusSlider.transform.Find("Text").GetComponent<Text>().text = "Particle radius = " + particleRadiusSlider.value.ToString("0.000");
        }
        if(vicsekLike)
        {
            gxSlider.transform.Find("Text").GetComponent<Text>().text = "Gx=" + gxSlider.value.ToString("0.000");
            gySlider.transform.Find("Text").GetComponent<Text>().text = "Gy=" + gySlider.value.ToString("0.000");
            viewRangeSlider.transform.Find("Text").GetComponent<Text>().text = "View range=" + viewRangeSlider.value.ToString("0.0");
            particleFovSlider.transform.Find("Text").GetComponent<Text>().text = "FOV=" + particleFovSlider.value.ToString("0.0");
            noiseAngleSlider.transform.Find("Text").GetComponent<Text>().text = "Noise angle=" + noiseAngleSlider.value.ToString("0.0");
            particlePropelSlider.transform.Find("Text").GetComponent<Text>().text = "Propulsion=" + particlePropelSlider.value.ToString("0.000");
        }
        if(applyChange.isOn)
        {
            if(vicsekLike)
            {
                gx = gxSlider.value;
                gy = gySlider.value;
                particleViewRange = viewRangeSlider.value;
                particleFov = particleFovSlider.value;
                noiseAngleDeg = noiseAngleSlider.value;
                particlePropulsion = particlePropelSlider.value;
            }
            if(isFaller)
            {
                gx = gxSlider.value;
                gy = gySlider.value;
                for (int i = 0; i < 3; i++)
                {
                    if(destModeCoeffs[i+1]!=CosSliders[i].value) destModeCoeffs[i+1] = CosSliders[i].value;
                    if(destModeSinCoeffs[i+1]!=SinSliders[i].value) destModeSinCoeffs[i+1] = SinSliders[i].value;
                }
            }
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
                                                                        + new Vector2((spawnPos[i + 0*newparticleCount]*plotSizeDelta.x)/DIM_X,(spawnPos[i + 1*newparticleCount]*plotSizeDelta.y)/DIM_Y))*plotScale.x;
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

        loopCountSlider.transform.Find("Text").GetComponent<Text>().text = "Steps per frame=" + ((int)loopCountSlider.value).ToString();
        for (int i = 0; i < (int)loopCountSlider.value; i++)
        {
            LBMStep();
        }
        UpdatePlot();
        if(isFaller) pTheta = polygonParticle.theta;
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
        if(!noParticles&&!isFaller)
        {
            for (int i = 0; i < particleCount; i++)
            {
                // roundParticles.PlotParticlePerimeter(ref plotPixels, particleColor);
                // roundParticles.PlotParticleFill(ref plotPixels);
                roundParticles.objs[i].GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
                                                                                    + new Vector2((roundParticles.pos[i + 0*particleCount]*plotSizeDelta.x)/DIM_X,(roundParticles.pos[i + 1*particleCount]*plotSizeDelta.y)/DIM_Y))*plotScale.x;
                roundParticles.objs[i].transform.rotation = Quaternion.Euler(0,0,((roundParticles.theta[i]*180f)/Mathf.PI));
            }
        }
        if(isFaller)
        {
            // polygonParticle.PlotParticlePerimeter(ref plotPixels,DIM_X,perimColor);
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
                fx[i + j*DIM_X] = forceFromParticlesX[i + j*DIM_X]; 
                if(isFaller||vicsekLike) fy[i + j*DIM_X] = forceFromParticlesY[i + j*DIM_X];
                else fy[i + j*DIM_X] = rbetag*(temp[i + j*DIM_X] - 0.5f) + forceFromParticlesY[i + j*DIM_X];
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i + j*DIM_X] + cy[k]*v[i + j*DIM_X];     
                    f0[k + (i + j*DIM_X)*9] = wf[k]*rho[i + j*DIM_X]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);

                    f[k + (i + j*DIM_X)*9] = f[k + (i + j*DIM_X)*9] - (f[k + (i + j*DIM_X)*9] - f0[k + (i + j*DIM_X)*9])/tauf;
                    f[k + (i + j*DIM_X)*9] += 3f*wf[k]*(cx[k]*fx[i + j*DIM_X] + cy[k]*fy[i + j*DIM_X]);
                    ftmp[k + (i + j*DIM_X)*9] = f[k + (i + j*DIM_X)*9];
                    if(k<5&&!isFaller&&!vicsekLike)
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
                        if(k<5&&!isFaller&&!vicsekLike)
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
            
            if(!isFaller&&!vicsekLike)
            {
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
            }
        }
        for(int i = 0; i < DIM_X; i++){
            f[4+ (i + (DIM_Y-1)*DIM_X)*9] = f[2+ (i + (DIM_Y-1)*DIM_X)*9];
            f[7+ (i + (DIM_Y-1)*DIM_X)*9] = f[5+ (i + (DIM_Y-1)*DIM_X)*9];
            f[8+ (i + (DIM_Y-1)*DIM_X)*9] = f[6+ (i + (DIM_Y-1)*DIM_X)*9]; 
            f[2+ (i + 0*DIM_X)*9] = f[4+ (i + 0*DIM_X)*9]; 
            f[5+ (i + 0*DIM_X)*9] = f[7+ (i + 0*DIM_X)*9]; 
            f[6+ (i + 0*DIM_X)*9] = f[8+ (i + 0*DIM_X)*9]; 

            if(!isFaller&&!vicsekLike)
            {
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
            
            if(vicsekLike&&align.isOn)
            {
                Vector2 averageMuki = new Vector2(0f,0f);
                int closeParticleCount = 0;
                for (int k = 0; k < particleCount; k++)
                {
                    if(k==n) continue;
                    for (int i = 0; i < 2; i++)
                    {
                        float dist = roundParticles.ParticleDistance(n,k);
                        if(dist < 2.0f*roundParticles.radius[n] + zeta)
                        {
                            roundParticles.forceFromCollisions[n + i*particleCount] += (roundParticles.pos[n + i*particleCount] - roundParticles.pos[k + i*particleCount])*(2.0f*roundParticles.radius[n] - dist + zeta)*(2.0f*roundParticles.radius[n] - dist + zeta)/epsw;
                        }
                        if(dist < 2.0f*roundParticles.radius[n] + particleViewRange)
                        {
                            averageMuki += new Vector2(Mathf.Sin(roundParticles.theta[k]),Mathf.Cos(roundParticles.theta[k]));
                            closeParticleCount++;
                        }
                    }
                }
                if(closeParticleCount!=0)
                {
                    averageMuki /= closeParticleCount;
                    roundParticles.theta[n] = Mathf.Atan2(averageMuki[0],averageMuki[1]);
                }
                roundParticles.theta[n] +=(Random.Range(-noiseAngleDeg,noiseAngleDeg)*Mathf.PI)/180f;
            }
            else
            {
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
                if(vicsekLike) roundParticles.theta[n] +=(Random.Range(-noiseAngleDeg,noiseAngleDeg)*Mathf.PI)/180f;
            }
            

            float wallDistance;
            float angle;
            if(vicsekLike)
            {
                wallDistance = Mathf.Abs(roundParticles.pos[n + 1*particleCount]);
                if(wallDistance < particleViewRange)
                {
                    if(wallDistance < roundParticles.radius[n]) roundParticles.pos[n + 1*particleCount] = roundParticles.radius[n];
                    angle = Vector2.SignedAngle(
                        new Vector2(0,-1),
                        new Vector2(-Mathf.Sin(roundParticles.theta[n]),Mathf.Cos(roundParticles.theta[n]))
                    );
                    if(Mathf.Abs(angle) < particleFov/2f)
                    {
                        float turnAngle = particleFov/2f - Mathf.Abs(angle);
                        if(angle < 0) turnAngle *= -1;
                        roundParticles.theta[n] += (turnAngle*Mathf.PI)/180f;
                    }
                }
                wallDistance = Mathf.Abs(DIM_Y-1-roundParticles.pos[n + 1*particleCount]);
                if(wallDistance < particleViewRange)
                {
                    if(wallDistance < roundParticles.radius[n]) roundParticles.pos[n + 1*particleCount] = DIM_Y-1-roundParticles.radius[n];
                    angle = Vector2.SignedAngle(
                        new Vector2(0,1),
                        new Vector2(-Mathf.Sin(roundParticles.theta[n]),Mathf.Cos(roundParticles.theta[n]))
                    );
                    if(Mathf.Abs(angle) < particleFov/2f)
                    {
                        float turnAngle = particleFov/2f - Mathf.Abs(angle);
                        if(angle < 0) turnAngle *= -1;
                        roundParticles.theta[n] +=(turnAngle*Mathf.PI)/180f;
                    }
                }
                wallDistance = Mathf.Abs(roundParticles.pos[n + 0*particleCount]);
                if(wallDistance < particleViewRange)
                {
                    if(wallDistance < roundParticles.radius[n]) roundParticles.pos[n + 0*particleCount] =roundParticles.radius[n];
                    angle = Vector2.SignedAngle(
                        new Vector2(-1,0),
                        new Vector2(-Mathf.Sin(roundParticles.theta[n]),Mathf.Cos(roundParticles.theta[n]))
                    );
                    if(Mathf.Abs(angle) < particleFov/2f)
                    {
                        float turnAngle = particleFov/2f - Mathf.Abs(angle);
                        if(angle < 0) turnAngle *= -1;
                        roundParticles.theta[n] += (turnAngle*Mathf.PI)/180f;
                    }
                }
                wallDistance = Mathf.Abs(DIM_X-1-roundParticles.pos[n + 0*particleCount]);
                if(wallDistance < particleViewRange)
                {
                    if(wallDistance < roundParticles.radius[n]) roundParticles.pos[n + 0*particleCount] = DIM_X-1-roundParticles.radius[n];
                    angle = Vector2.SignedAngle(
                        new Vector2(1,0),
                        new Vector2(-Mathf.Sin(roundParticles.theta[n]),Mathf.Cos(roundParticles.theta[n]))
                    );
                    if(Mathf.Abs(angle) < particleFov/2f)
                    {
                        float turnAngle = particleFov/2f - Mathf.Abs(angle);
                        if(angle < 0) turnAngle *= -1;
                        roundParticles.theta[n] +=(turnAngle*Mathf.PI)/180f;
                    }
                }
            }
            else
            {
                wallDistance = Mathf.Abs(roundParticles.pos[n + 1*particleCount]);
                if(wallDistance < roundParticles.radius[n]) roundParticles.pos[n + 1*particleCount] = roundParticles.radius[n];
                wallDistance = Mathf.Abs(DIM_Y-1-roundParticles.pos[n + 1*particleCount]);
                if(wallDistance < roundParticles.radius[n]) roundParticles.pos[n + 1*particleCount] = DIM_Y-1-roundParticles.radius[n];
                wallDistance = Mathf.Abs(roundParticles.pos[n + 0*particleCount]);
                if(wallDistance < roundParticles.radius[n]) roundParticles.pos[n + 0*particleCount] =roundParticles.radius[n];
                wallDistance = Mathf.Abs(DIM_X-1-roundParticles.pos[n + 0*particleCount]);
                if(wallDistance < roundParticles.radius[n]) roundParticles.pos[n + 0*particleCount] = DIM_X-1-roundParticles.radius[n];
            }
            roundParticles.forceFromFluid[n + 0*particleCount] = 0f;
            roundParticles.forceFromFluid[n + 1*particleCount] = 0f;
            roundParticles.torque[n] = 0f;
            for(int m = 0; m < roundParticles.perimeterPointCount[n] ; m++) 
            {
                roundParticles.perimeterFluidVel[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] = 0f;
                roundParticles.perimeterFluidVel[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] = 0f;
                angle = Vector2.Angle(
                    new Vector2(roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] - roundParticles.pos[n + 0*particleCount],roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] - roundParticles.pos[n + 1*particleCount]),
                    new Vector2(-Mathf.Sin(roundParticles.theta[n]),Mathf.Cos(roundParticles.theta[n]))
                );
                if(!particleConstantNoTemp)
                {
                    
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

                if(angle>90f&&vicsekLike)
                {
                    Vector2 vect2Center = new Vector2(
                        roundParticles.pos[n + 0*particleCount] - roundParticles.perimeterPos[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount],
                        roundParticles.pos[n + 1*particleCount] - roundParticles.perimeterPos[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount]
                    );
                    vect2Center.Normalize();
                    roundParticles.forceOnPerimeter[n + (m + 0*roundParticles.maxPerimeterPointCount)*particleCount] += particlePropulsion * vect2Center[0] * Mathf.Cos((angle*Mathf.PI)/180f);
                    roundParticles.forceOnPerimeter[n + (m + 1*roundParticles.maxPerimeterPointCount)*particleCount] += particlePropulsion * vect2Center[1] * Mathf.Cos((angle*Mathf.PI)/180f);
                }
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
            
            if(vicsekLike)roundParticles.UpdatePosVel(n,new float[2]{gx,gy});
            else roundParticles.UpdatePosVel(n,gravity);
            roundParticles.UpdateOmegaTheta(n);
            roundParticles.UpdatePerimeter(n);
        }
    }

    void FallerImmersedBoundary()
    {
        float tmp1,tmp2,tmp3;
        polygonParticle.forceFromFluid[0] = 0f;
        polygonParticle.forceFromFluid[1] = 0f;
        polygonParticle.torque = 0f;

        float N,S,W,E;
        N = 0;
        S = Mathf.Infinity;
        E = 0;
        W = Mathf.Infinity;
        for(int m = 0; m < polygonParticle.perimeterPointCount ; m++) 
        {
            N = Mathf.Max(N,polygonParticle.perimeterPos[m,1]);
            S = Mathf.Min(S,polygonParticle.perimeterPos[m,1]);
            E = Mathf.Max(E,polygonParticle.perimeterPos[m,0]);
            W = Mathf.Min(W,polygonParticle.perimeterPos[m,0]);
        }

        if(N>DIM_Y) polygonParticle.pos[1] -= N - DIM_Y;
        if(S<0)     polygonParticle.pos[1] += -S;
        if(E>DIM_X) polygonParticle.pos[0] -= E - DIM_X;
        if(W<0)     polygonParticle.pos[0] += -W;

        // 固体表面の流体の速度を計算
        for(int m = 0; m < polygonParticle.perimeterPointCount ; m++) 
        {
            // uet[n,m] = 0.0f; vet[n,m] = 0.0f;
            polygonParticle.perimeterFluidVel[m,0] = 0f;
            polygonParticle.perimeterFluidVel[m,1] = 0f;
            for(int i = (int)polygonParticle.perimeterPos[m,0] - 3; i < (int)polygonParticle.perimeterPos[m,0] + 3; i++)
            {
                for(int j = (int)polygonParticle.perimeterPos[m,1] - 3; j < (int)polygonParticle.perimeterPos[m,1] + 3; j++)
                {
                    tmp1 = Mathf.Abs(polygonParticle.perimeterPos[m,0] - (float)i);
                    tmp2 = Mathf.Abs(polygonParticle.perimeterPos[m,1] - (float)j);
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
                    if((j<DIM_Y&&j>=0)&&(i<DIM_X&&i>=0))
                    {
                        polygonParticle.perimeterFluidVel[m,0] += u[i + j*DIM_X]*tmp3;
                        polygonParticle.perimeterFluidVel[m,1] += v[i + j*DIM_X]*tmp3;
                    }
                } 
            }
            polygonParticle.forceOnPerimeter[m,0] = polygonParticle.perimeterVel[m,0] - polygonParticle.perimeterFluidVel[m,0];
            polygonParticle.forceOnPerimeter[m,1] = polygonParticle.perimeterVel[m,1] - polygonParticle.perimeterFluidVel[m,1];

            // 固体が外部に与える力を計算
            for(int i = (int)polygonParticle.perimeterPos[m,0] - 3; i < (int)polygonParticle.perimeterPos[m,0] + 3; i++)
            {
                for(int j = (int)polygonParticle.perimeterPos[m,1] - 3; j < (int)polygonParticle.perimeterPos[m,1] + 3; j++)
                {
                    tmp1 = Mathf.Abs(polygonParticle.perimeterPos[m,0] - (float)i);
                    tmp2 = Mathf.Abs(polygonParticle.perimeterPos[m,1] - (float)j);
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
                    if((j<DIM_Y&&j>=0)&&(i<DIM_X&&i>=0))
                    {
                        forceFromParticlesX[i + j*DIM_X] += polygonParticle.forceOnPerimeter[m,0] * tmp3 * 0.5f;
                        forceFromParticlesY[i + j*DIM_X] += polygonParticle.forceOnPerimeter[m,1] * tmp3 * 0.5f;
                    }
                } 
            }
            polygonParticle.forceFromFluid[0] += polygonParticle.forceOnPerimeter[m,0];
            polygonParticle.forceFromFluid[1] += polygonParticle.forceOnPerimeter[m,1];
            polygonParticle.torque += polygonParticle.forceOnPerimeter[m,1] * (polygonParticle.perimeterPos[m,0] - polygonParticle.pos[0]) 
                                    - polygonParticle.forceOnPerimeter[m,0] * (polygonParticle.perimeterPos[m,1] - polygonParticle.pos[1]);
        } 

        polygonParticle.forceFromFluid[0] *= -0.5f;  
        polygonParticle.forceFromFluid[1] *= -0.5f; 
        polygonParticle.torque *= -0.5f; 

        polygonParticle.UpdatePosVel(new float[2]{gx,gy});
        polygonParticle.UpdateOmegaTheta();
        UpdatePolygonPerimeter();
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
        tracePoints.Clear();
        DIM_X = int.Parse(inputFieldW.text);
        DIM_Y = int.Parse(inputFieldH.text);
        time = 0;
        plotTexture = new Texture2D(DIM_X,DIM_Y);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),UnityEngine.Vector2.zero);

        if(DIM_X>DIM_Y)
        {
            plotImage.transform.GetComponent<RectTransform>().sizeDelta = new Vector2(
                originalPlotSizeDelta.x,
                (originalPlotSizeDelta.x*(float)DIM_Y)/(float)DIM_X
            );
        }
        else if(DIM_X<DIM_Y)
        {
            plotImage.transform.GetComponent<RectTransform>().sizeDelta = new Vector2(
                (originalPlotSizeDelta.y*(float)DIM_X)/(float)DIM_Y,
                originalPlotSizeDelta.y
            );
        }
        else
        {
            plotImage.transform.GetComponent<RectTransform>().sizeDelta = originalPlotSizeDelta;
        }
        vectorField.DestroyAll();
        vectorField.aspectRatio = (float)DIM_X/(float)DIM_Y;
        vectorField.Start();
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
        plotSizeDelta = plotImage.transform.GetComponent<RectTransform>().sizeDelta;
        plotScale = plotImage.transform.localScale;
        if(!isFaller)
        {
            float[] spawnPos = new float[particleCount*2];
            GameObject[] objs = new GameObject[particleCount];
            particleRadiusSlider.value = particleRadius;
            float scale = (2f*particleRadius*plotImage.transform.GetComponent<RectTransform>().sizeDelta.x * plotImage.transform.localScale.x)/(DIM_X*particlePrefab.GetComponent<RectTransform>().sizeDelta.x);
            particleCountSlider.value = particleCount;
            for (int i = 0; i < particleCount; i++)
            {
                Destroy(roundParticles.objs[i]);
                spawnPos[i + 0*particleCount] = UnityEngine.Random.Range(particleRadius,DIM_X-particleRadius);
                spawnPos[i + 1*particleCount] = UnityEngine.Random.Range(particleRadius,DIM_Y-particleRadius);
                objs[i] = Instantiate(particlePrefab,particleParent);
                objs[i].GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
                                                                + new Vector2((spawnPos[i + 0*particleCount]*plotSizeDelta.x)/DIM_X,(spawnPos[i + 1*particleCount]*plotSizeDelta.y)/DIM_Y))*plotScale.x;
                objs[i].transform.localScale = new Vector3(scale,scale,1);
            }
            roundParticles = new RoundParticlesForCompute(DIM_X,particleCount,particleDensity,particleRadius,spawnPos,objs,Random.Range(0f,2f*Mathf.PI));
        }
        else
        {
            UpdatePolygonPerimeter(true);
        }

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

    void UpdatePolygonPerimeter(bool initParticle = false)
    {
        List<Vector3> points = new List<Vector3>();
        float circumferenceProgressPerStep = (float)1/lineRes;
        float TAU = 2*Mathf.PI;
        float radianProgressPerStep = circumferenceProgressPerStep*TAU;
        float scale;
        float sqrdSum = 0f;
        float polygonTheta = 0f;
        area = sizeSlider.value;
        sizeSlider.transform.Find("Text").GetComponent<Text>().text = "Size = " + sizeSlider.value.ToString("0.0");
        if(initParticle) polygonTheta = 0f;
        else
        {
            polygonTheta = polygonParticle.theta;
            if(area != polygonParticle.volume) polygonParticle.UpdateVolume(area);
        }
        for (int i = 0; i < lineRes; i++)
        {
            float currentRadian = radianProgressPerStep*i + polygonTheta;
            float newRadius = 0f;
            for (int j = 0; j < modes; j++)
            {
                newRadius += modeCoeffs[j]*Mathf.Cos(j*(currentRadian)) + modeSinCoeffs[j]*Mathf.Sin(j*(currentRadian));
            }
            sqrdSum += newRadius * newRadius;
        }
        scale = Mathf.Sqrt((2*area)/(sqrdSum * lineRes * Mathf.Sin(TAU/lineRes)));

        float perimeterLength = 2f * Mathf.PI * modeCoeffs[0] * scale;
        float segmentLength = 1f/(segmentDivisor* (DIM_X/((canvas.localScale.x)*(plotImage.transform.GetComponent<RectTransform>().sizeDelta.x*plotImage.transform.localScale.x))));
        float theta = 0f;
        float dtheta,r;
        int spawnedDots = 0;
        Vector3 spawnPos;
        List<float> thetaList = new List<float>();
        while (theta <= 2*Mathf.PI)
        {
            r = 0f;
            for (int i = 0; i < modes; i++)
            {
                r += modeCoeffs[i]*Mathf.Cos(i*(theta)) + modeSinCoeffs[i]*Mathf.Sin(i*(theta));
            }
            r *= scale;
            if(r <= Mathf.Epsilon)
            {
                theta += 0.1f * Mathf.PI;
                continue;
            }
            spawnPos = new Vector3(Mathf.Cos((theta+polygonTheta)), Mathf.Sin((theta+polygonTheta)),0);
            spawnPos *= r;
            thetaList.Add(theta);
            if(dotPositions.Count < spawnedDots + 1)
            {
                dotPositions.Add(spawnPos);
            }
            else
            {
                dotPositions[spawnedDots] = spawnPos;
            }
            spawnedDots++;

            dtheta = segmentLength/r;
            theta += dtheta;
        }
        if(spawnedDots<dotPositions.Count)
        {
            int a = dotPositions.Count;
            for (int i = 0; i < a - spawnedDots; i++)
            {
                dotPositions.RemoveAt(dotPositions.Count-1);
            }
        }
        if(initParticle){
            polygonParticle = new PolygonParticle(spawnedDots,particleDensity,area,new float[2]{(float)(DIM_X)/2f,((float)(DIM_Y)/2f)});
        }
        polygonParticle.perimeterPointCount =  dotPositions.Count;
        polygonParticle.perimeterPos = new float[polygonParticle.perimeterPointCount,2];
        List<Vector3> linePoints = new List<Vector3>();
        polygonParticle.forceOnPerimeter = new float[polygonParticle.perimeterPointCount,2];
        polygonParticle.perimeterFluidVel = new float[polygonParticle.perimeterPointCount,2];
        polygonParticle.perimeterVel = new float[polygonParticle.perimeterPointCount,2];
        List<int> triangleList = new List<int>();
        for (int i = 0; i < polygonParticle.perimeterPointCount; i++)
        {
            if(i<polygonParticle.perimeterPointCount-1)
            {
                triangleList.Add(polygonParticle.perimeterPointCount);
                triangleList.Add(i);
                triangleList.Add(i+1);
            }
            else
            {
                triangleList.Add(polygonParticle.perimeterPointCount);
                triangleList.Add(i);
                triangleList.Add(0);
            }
            float newposx = dotPositions[i].x/canvas.localScale.x;
            newposx *= DIM_X/(plotImage.transform.GetComponent<RectTransform>().sizeDelta.x*plotImage.transform.localScale.x);
            float newposy = dotPositions[i].y/canvas.localScale.y;
            newposy *= DIM_Y/(plotImage.transform.GetComponent<RectTransform>().sizeDelta.y*plotImage.transform.localScale.y);
            polygonParticle.perimeterPos[i,0] = polygonParticle.pos[0] + newposx;
            polygonParticle.perimeterPos[i,1] = polygonParticle.pos[1] + newposy;
            linePoints.Add(new Vector3(
                (((polygonParticle.perimeterPos[i,0] - (float)DIM_X/2f)/((float)DIM_X/2f))*plotSizeDelta.x*plotScale.x/2f + plotAnchPos.x)*canvas.localScale.x,
                (((polygonParticle.perimeterPos[i,1] - (float)DIM_Y/2f)/((float)DIM_Y/2f))*plotSizeDelta.y*plotScale.y/2f + plotAnchPos.y)*canvas.localScale.y,
                0
            ));
            polygonParticle.perimeterVel[i,0] = polygonParticle.vel[0] - polygonParticle.omega*(polygonParticle.perimeterPos[i,1] - polygonParticle.pos[1]);
            polygonParticle.perimeterVel[i,1] = polygonParticle.vel[1] + polygonParticle.omega*(polygonParticle.perimeterPos[i,0] - polygonParticle.pos[0]);
            
            Vector2 newPosVec = new Vector2(newposx, newposy);
            newPosVec.Normalize();
            for (int j = 1; j < 4; j++)
            {
                float cosCoeffVel = 0f;
                float sinCoeffVel = 0f;
                float cosDiff,sinDiff;
                cosDiff = originalModeCoeffs[j]-destModeCoeffs[j];
                sinDiff = originalModeSinCoeffs[j]-destModeSinCoeffs[j];
                if(CosShuukiSliders[j-1].transform.Find("Toggle").GetComponent<Toggle>().isOn)
                {
                    cosCoeffVel += (-originalModeCoeffs[j]*(2f*Mathf.PI/CosShuukiSliders[j-1].value))*Mathf.Sin(2f*Mathf.PI*((float)time/CosShuukiSliders[j-1].value));
                    if(cosDiff!=0f) cosCoeffVel +=(-cosDiff/Mathf.Abs(cosDiff))*coeffChangeSpeed*Mathf.Cos(2f*Mathf.PI*((float)time/CosShuukiSliders[j-1].value));
                    polygonParticle.perimeterVel[i,0] += scale*newPosVec[0] * cosCoeffVel * Mathf.Cos(j*(thetaList[i]));
                    polygonParticle.perimeterVel[i,1] += scale*newPosVec[1] * cosCoeffVel * Mathf.Cos(j*(thetaList[i]));
                }
                else
                {
                    if(cosDiff!=0f)cosCoeffVel += (-cosDiff/Mathf.Abs(cosDiff))*coeffChangeSpeed;
                    polygonParticle.perimeterVel[i,0] += scale*newPosVec[0] * cosCoeffVel * Mathf.Cos(j*(thetaList[i]));
                    polygonParticle.perimeterVel[i,1] += scale*newPosVec[1] * cosCoeffVel * Mathf.Cos(j*(thetaList[i]));
                }
                if(SinShuukiSliders[j-1].transform.Find("Toggle").GetComponent<Toggle>().isOn)
                {
                    sinCoeffVel += -originalModeSinCoeffs[j]*(2f*Mathf.PI/SinShuukiSliders[j-1].value)*Mathf.Sin(2f*Mathf.PI*((float)time/SinShuukiSliders[j-1].value));
                    if(sinDiff!=0f) sinCoeffVel += (-sinDiff/Mathf.Abs(sinDiff))*coeffChangeSpeed*Mathf.Cos(2f*Mathf.PI*((float)time/SinShuukiSliders[j-1].value));
                    polygonParticle.perimeterVel[i,0] += scale*newPosVec[0] * sinCoeffVel * Mathf.Sin(j*(thetaList[i]));
                    polygonParticle.perimeterVel[i,1] += scale*newPosVec[1] * sinCoeffVel * Mathf.Sin(j*(thetaList[i]));
                }
                else
                {
                    if(sinDiff!=0f)sinCoeffVel += (-sinDiff/Mathf.Abs(sinDiff))*coeffChangeSpeed;
                    polygonParticle.perimeterVel[i,0] += scale*newPosVec[0] * sinCoeffVel * Mathf.Sin(j*(thetaList[i]));
                    polygonParticle.perimeterVel[i,1] += scale*newPosVec[1] * sinCoeffVel * Mathf.Sin(j*(thetaList[i]));
                }
            }
        }

        // line.positionCount = linePoints.Count;
        // line.SetPositions(linePoints.ToArray());
        mesh.Clear();
        linePoints.Add(new Vector3(
            (((polygonParticle.pos[0] - (float)DIM_X/2f)/((float)DIM_X/2f))*plotSizeDelta.x*plotScale.x/2f + plotAnchPos.x)*canvas.localScale.x,
            (((polygonParticle.pos[1] - (float)DIM_Y/2f)/((float)DIM_Y/2f))*plotSizeDelta.y*plotScale.y/2f + plotAnchPos.y)*canvas.localScale.y,
            0
        ));
        mesh.vertices = linePoints.ToArray();
        mesh.triangles = triangleList.ToArray();
        if(tracePoints.Count == 0)
        {
            tracePoints.Add(new Vector3(
                (((polygonParticle.pos[0] - (float)DIM_X/2f)/((float)DIM_X/2f))*plotSizeDelta.x*plotScale.x/2f + plotAnchPos.x)*canvas.localScale.x,
                (((polygonParticle.pos[1] - (float)DIM_Y/2f)/((float)DIM_Y/2f))*plotSizeDelta.y*plotScale.y/2f + plotAnchPos.y)*canvas.localScale.y,
                0
            ));
            tracePoints.Add(new Vector3(
                (((polygonParticle.pos[0] - (float)DIM_X/2f)/((float)DIM_X/2f))*plotSizeDelta.x*plotScale.x/2f + plotAnchPos.x)*canvas.localScale.x,
                (((polygonParticle.pos[1] - (float)DIM_Y/2f)/((float)DIM_Y/2f))*plotSizeDelta.y*plotScale.y/2f + plotAnchPos.y)*canvas.localScale.y,
                0
            ));
            traceLine.positionCount = tracePoints.Count;
            traceLine.SetPositions(tracePoints.ToArray());
        }
        else if((tracePoints[tracePoints.Count-1] - tracePoints[tracePoints.Count-2]).sqrMagnitude > traceThres*traceThres)
        {
            tracePoints.Add(new Vector3(
                (((polygonParticle.pos[0] - (float)DIM_X/2f)/((float)DIM_X/2f))*plotSizeDelta.x*plotScale.x/2f + plotAnchPos.x)*canvas.localScale.x,
                (((polygonParticle.pos[1] - (float)DIM_Y/2f)/((float)DIM_Y/2f))*plotSizeDelta.y*plotScale.y/2f + plotAnchPos.y)*canvas.localScale.y,
                0
            ));
            traceLine.positionCount = tracePoints.Count;
            traceLine.SetPositions(tracePoints.ToArray());
        }
        else
        {
            tracePoints[tracePoints.Count-1] = new Vector3(
                (((polygonParticle.pos[0] - (float)DIM_X/2f)/((float)DIM_X/2f))*plotSizeDelta.x*plotScale.x/2f + plotAnchPos.x)*canvas.localScale.x,
                (((polygonParticle.pos[1] - (float)DIM_Y/2f)/((float)DIM_Y/2f))*plotSizeDelta.y*plotScale.y/2f + plotAnchPos.y)*canvas.localScale.y,
                0
            );
            traceLine.SetPositions(tracePoints.ToArray());
        }
    }

    public void OnShindoAllChange()
    {
        for (int i = 0; i < 3; i++)
        {
            CosShuukiSliders[i].transform.Find("Toggle").GetComponent<Toggle>().isOn = shindoAll.isOn;
            SinShuukiSliders[i].transform.Find("Toggle").GetComponent<Toggle>().isOn = shindoAll.isOn;
        }
    }

    public void AllZero()
    {
        for (int i = 0; i < 3; i++)
        {
            CosSliders[i].value = 0f;
            SinSliders[i].value = 0f;
        }
        gxSlider.value = 0f;
        gySlider.value = 0f;
    }

    public void ZeroGrav()
    {
        gxSlider.value = 0f;
        gySlider.value = 0f;
    }

    public void BackToMenu()
    {
        SceneManager.LoadScene("Menu");
    }
}