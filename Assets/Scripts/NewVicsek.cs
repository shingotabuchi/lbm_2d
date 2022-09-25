using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class NewVicsek : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public VectorField vectorField;
    public HeatMapMode mode = HeatMapMode.Speed;
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
    public float tauf = 0.8f;
    public int particleCount = 25;
    public float particleDensity = 1.25f;
    public float particleRadius = 1.25f;
    public float zeta = 1f;
    public float epsw = 1f;
    RoundParticle[] roundParticles;

    public GameObject particlePrefab;
    public Transform particleParent;
    Vector2 plotAnchPos; 
    Vector2 plotSizeDelta; 
    Vector3 plotScale; 

    public Toggle applyChange;
    public Toggle normalize;

    public GameObject heatMap;

    Vector2 heatMapOriginalPos;
    float heatMapInitHeight;
    float heatMapHeight;

    public Transform canvas;

    public Slider particleRadiusSlider;

    public Slider particleCountSlider;

    public float gx = 0f;
    public float gy = 0f;

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
    public float particleViewRange,particleFov,noiseAngleDeg,particlePropulsion;
    // Start is called before the first frame update
    void Start()
    {
        particleViewRange = viewRangeSlider.value;
        particleFov = particleFovSlider.value;
        noiseAngleDeg = noiseAngleSlider.value;
        particlePropulsion = particlePropelSlider.value;
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
        taug = 3.0f*chi + 0.5f;

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
        plotAnchPos = plotImage.transform.GetComponent<RectTransform>().anchoredPosition;
        plotSizeDelta = plotImage.transform.GetComponent<RectTransform>().sizeDelta;
        plotScale = plotImage.transform.localScale;
        roundParticles = new RoundParticle[particleCount];
        float scale = (2f*particleRadius*plotImage.transform.GetComponent<RectTransform>().sizeDelta.x * plotImage.transform.localScale.x)/(DIM_X*particlePrefab.GetComponent<RectTransform>().sizeDelta.x);
        for (int i = 0; i < particleCount; i++)
        {
            float[] spawnPos =  new float[2]{UnityEngine.Random.Range(particleRadius,DIM_X-particleRadius),UnityEngine.Random.Range(particleRadius,DIM_Y-particleRadius)};
            roundParticles[i] = new RoundParticle(particleDensity,particleRadius,spawnPos,Random.Range(0f,2f*Mathf.PI));
            roundParticles[i].obj = Instantiate(particlePrefab,particleParent);
            roundParticles[i].obj.GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
                                                                                + new Vector2((spawnPos[0]*plotSizeDelta.x)/DIM_X,(spawnPos[1]*plotSizeDelta.y)/DIM_Y))*plotScale.x;
            roundParticles[i].obj.transform.localScale = new Vector3(scale,scale,1);
        }
        
        maxSpeed = 0f;
        minSpeed = Mathf.Infinity;
        maxRho = 0f;
        minRho = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u[i + j*DIM_X] = 0.0f; v[i + j*DIM_X] = 0.0f; 
                rho[i + j*DIM_X] = 1.0f;

                u2 = u[i + j*DIM_X]*u[i + j*DIM_X] + v[i + j*DIM_X]*v[i + j*DIM_X];   

                forceFromParticlesX[i + j*DIM_X] = 0f;
                forceFromParticlesY[i + j*DIM_X] = 0f;
                
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i + j*DIM_X] + cy[k]*v[i + j*DIM_X];     
                    f0[k + (i + j*DIM_X)*9] = wf[k]*rho[i + j*DIM_X]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
                    f[k + (i + j*DIM_X)*9] = f0[k + (i + j*DIM_X)*9];
                }
                speed[i + j*DIM_X] = Mathf.Sqrt(u2);
                maxSpeed = Mathf.Max(maxSpeed,speed[i + j*DIM_X]);
                minSpeed = Mathf.Min(minSpeed,speed[i + j*DIM_X]);
                maxRho = Mathf.Max(maxRho,rho[i + j*DIM_X]);
                minRho = Mathf.Min(minRho,rho[i + j*DIM_X]);
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
        ImmersedBoundary();
    }

    // Update is called once per frame
    void Update()
    {
        normalizeHeatMap = normalize.isOn;
        mode = HeatMapMode.Speed;

        gxSlider.transform.Find("Text").GetComponent<Text>().text = "Gx=" + gxSlider.value.ToString("0.000");
        gySlider.transform.Find("Text").GetComponent<Text>().text = "Gy=" + gySlider.value.ToString("0.000");
        viewRangeSlider.transform.Find("Text").GetComponent<Text>().text = "View range=" + viewRangeSlider.value.ToString("0.0");
        particleFovSlider.transform.Find("Text").GetComponent<Text>().text = "FOV=" + particleFovSlider.value.ToString("0.0");
        noiseAngleSlider.transform.Find("Text").GetComponent<Text>().text = "Noise angle=" + noiseAngleSlider.value.ToString("0.0");
        particlePropelSlider.transform.Find("Text").GetComponent<Text>().text = "Propulsion=" + particlePropelSlider.value.ToString("0.000");

        if(applyChange.isOn)
        {
            gx = gxSlider.value;
            gy = gySlider.value;
            particleViewRange = viewRangeSlider.value;
            particleFov = particleFovSlider.value;
            noiseAngleDeg = noiseAngleSlider.value;
            particlePropulsion = particlePropelSlider.value;

            // if((int)particleCountSlider.value != particleCount)
            // {
            //     int newparticleCount = (int)particleCountSlider.value;
            //     float[] spawnPos = new float[newparticleCount*2];
            //     float[] newTheta = new float[newparticleCount];
            //     GameObject[] objs = new GameObject[newparticleCount];
            //     float scale = (2f*particleRadius*plotImage.transform.GetComponent<RectTransform>().sizeDelta.x * plotImage.transform.localScale.x)/(DIM_X*particlePrefab.GetComponent<RectTransform>().sizeDelta.x);
            //     for (int i = 0; i < newparticleCount; i++)
            //     {
                    
            //         if(i<particleCount)
            //         {
            //             newTheta[i] = roundParticles.theta[i];
            //             spawnPos[i + 0*newparticleCount] = roundParticles.pos[i + 0*particleCount];
            //             spawnPos[i + 1*newparticleCount] = roundParticles.pos[i + 1*particleCount];
            //             objs[i] = roundParticles.objs[i];
            //         }
            //         else
            //         {
            //             newTheta[i] = Random.Range(0f,2f*Mathf.PI);
            //             spawnPos[i + 0*newparticleCount] = UnityEngine.Random.Range(particleRadius,DIM_X-particleRadius);
            //             spawnPos[i + 1*newparticleCount] = UnityEngine.Random.Range(particleRadius,DIM_Y-particleRadius);
            //             objs[i] = Instantiate(particlePrefab,particleParent);
            //             objs[i].GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
            //                                                             + new Vector2((spawnPos[i + 0*newparticleCount]*plotSizeDelta.x)/DIM_X,(spawnPos[i + 1*newparticleCount]*plotSizeDelta.y)/DIM_Y))*plotScale.x;
            //             objs[i].transform.localScale = new Vector3(scale,scale,1);
            //         }
            //     }
            //     if(newparticleCount < particleCount)
            //     {
            //         for (int i = newparticleCount; i < particleCount; i++)
            //         {
            //             Destroy(roundParticles.objs[i]);
            //         }
            //     }
            //     particleCount = newparticleCount;
            //     roundParticles.UpdateParticleCount(particleCount,spawnPos,newTheta,objs);
            // }
            // if(particleRadiusSlider.value != particleRadius)
            // {
            //     particleRadius = particleRadiusSlider.value;
            //     float scale = (2f*particleRadius*plotImage.transform.GetComponent<RectTransform>().sizeDelta.x * plotImage.transform.localScale.x)/(DIM_X*particlePrefab.GetComponent<RectTransform>().sizeDelta.x);
            //     for (int i = 0; i < particleCount; i++)
            //     {
            //         roundParticles.objs[i].transform.localScale = new Vector3(scale,scale,1);
            //     }
            //     roundParticles.UpdateRadius(particleRadius);
            // }
        }

        loopCountSlider.transform.Find("Text").GetComponent<Text>().text = "Steps per frame=" + ((int)loopCountSlider.value).ToString();
        for (int i = 0; i < (int)loopCountSlider.value; i++)
        {
            LBMStep();
        }
        UpdatePlot();
    }

    void UpdatePlot()
    {
        if(maxSpeed == 0f)maxSpeed = 1f;
        if(maxTemp == 0f)maxTemp = 1f;
        // UpdateHeatMap();
        for (int i = 0; i < plotPixels.Length; i++)
        {
            plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X+(i/DIM_X)*DIM_X]-minSpeed,maxSpeed-minSpeed);
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
        for (int i = 0; i < particleCount; i++)
        {
            roundParticles[i].obj.GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
                                                                                + new Vector2((roundParticles[i].pos[0]*plotSizeDelta.x)/DIM_X,(roundParticles[i].pos[1]*plotSizeDelta.y)/DIM_Y))*plotScale;
            roundParticles[i].obj.transform.rotation = Quaternion.Euler(0,0,(roundParticles[i].theta*180f)/Mathf.PI);
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
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i + j*DIM_X] + cy[k]*v[i + j*DIM_X];     
                    f0[k + (i + j*DIM_X)*9] = wf[k]*rho[i + j*DIM_X]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);

                    f[k + (i + j*DIM_X)*9] = f[k + (i + j*DIM_X)*9] - (f[k + (i + j*DIM_X)*9] - f0[k + (i + j*DIM_X)*9])/tauf;
                    f[k + (i + j*DIM_X)*9] += 3f*wf[k]*(cx[k]*forceFromParticlesX[i + j*DIM_X] + cy[k]*forceFromParticlesY[i + j*DIM_X]);
                    // ftmp[k + (i + j*DIM_X)*9] = f[k + (i + j*DIM_X)*9];
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
        // ftmp = (float[])(f.Clone());
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
                    }
                } 
            }
        }
    }

    void Boundaries()
    {
        for (int j = 0; j < DIM_Y; j++)
        {
            f[1 + (0 + j*DIM_X)*9] = f[3 + (0 + j*DIM_X)*9];
            f[5 + (0 + j*DIM_X)*9] = f[7 + (0 + j*DIM_X)*9];
            f[8 + (0 + j*DIM_X)*9] = f[6 + (0 + j*DIM_X)*9];
            f[3 + (DIM_X-1 + j*DIM_X)*9] = f[1 + (DIM_X-1 + j*DIM_X)*9]; 
            f[7 + (DIM_X-1 + j*DIM_X)*9] = f[5 + (DIM_X-1 + j*DIM_X)*9]; 
            f[6 + (DIM_X-1 + j*DIM_X)*9] = f[8 + (DIM_X-1 + j*DIM_X)*9]; 
        }
        for(int i = 0; i < DIM_X; i++){
            f[4+ (i + (DIM_Y-1)*DIM_X)*9] = f[2+ (i + (DIM_Y-1)*DIM_X)*9];
            f[7+ (i + (DIM_Y-1)*DIM_X)*9] = f[5+ (i + (DIM_Y-1)*DIM_X)*9];
            f[8+ (i + (DIM_Y-1)*DIM_X)*9] = f[6+ (i + (DIM_Y-1)*DIM_X)*9]; 
            f[2+ (i + 0*DIM_X)*9] = f[4+ (i + 0*DIM_X)*9]; 
            f[5+ (i + 0*DIM_X)*9] = f[7+ (i + 0*DIM_X)*9]; 
            f[6+ (i + 0*DIM_X)*9] = f[8+ (i + 0*DIM_X)*9]; 
        }
    }

    void UpdateSpeedAndTemperature()
    {
        maxSpeed = 0f;
        minSpeed = Mathf.Infinity;
        maxRho = 0f;
        minRho = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u[i + j*DIM_X] = 0f; v[i + j*DIM_X] = 0f;
                rho[i + j*DIM_X] = f[0+(i+j*DIM_X)*9]; 
                for(int k = 1; k <= 8; k++)
                {
                    rho[i + j*DIM_X] = rho[i + j*DIM_X] + f[k + (i + j*DIM_X)*9];
                    u[i + j*DIM_X] =   u[i + j*DIM_X] + f[k + (i + j*DIM_X)*9]*cx[k];
                    v[i + j*DIM_X] =   v[i + j*DIM_X] + f[k + (i + j*DIM_X)*9]*cy[k];
                } 
                u[i + j*DIM_X] = u[i + j*DIM_X]/rho[i + j*DIM_X];
                v[i + j*DIM_X] = v[i + j*DIM_X]/rho[i + j*DIM_X];
                speed[i + j*DIM_X] = Mathf.Sqrt(u[i + j*DIM_X]*u[i + j*DIM_X] + v[i + j*DIM_X]*v[i + j*DIM_X]);
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
            roundParticles[n].forceFromCollisions[0] = 0f;
            roundParticles[n].forceFromCollisions[1] = 0f;
            // tmp1 = Mathf.Abs(roundParticles[n].pos[1] + roundParticles[n].radius); 
            // if(tmp1 < 2.0f*roundParticles[n].radius + zeta){
            //     roundParticles[n].forceFromCollisions[1] = (roundParticles[n].pos[1] + roundParticles[n].radius)*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
            // }
            // tmp1 = Mathf.Abs(DIM_Y-1-roundParticles[n].pos[1] + roundParticles[n].radius); 
            // if(tmp1 < 2.0f*roundParticles[n].radius + zeta){
            //     roundParticles[n].forceFromCollisions[1] = -(DIM_Y-1-roundParticles[n].pos[1] + roundParticles[n].radius)*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
            // }
            // tmp1 = Mathf.Abs(roundParticles[n].pos[0] + roundParticles[n].radius); 
            // if(tmp1 < 2.0f*roundParticles[n].radius + zeta){
            //     roundParticles[n].forceFromCollisions[0] = (roundParticles[n].pos[0] + roundParticles[n].radius)*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
            // }
            // tmp1 = Mathf.Abs(DIM_X-1-roundParticles[n].pos[0] + roundParticles[n].radius); 
            // if(tmp1 < 2.0f*roundParticles[n].radius + zeta){
            //     roundParticles[n].forceFromCollisions[0] = -(DIM_X-1-roundParticles[n].pos[0] + roundParticles[n].radius)*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
            // }

            if(align){
                Vector2 averageMuki = new Vector2(0f,0f);
                int closeParticleCount = 0;
                for (int k = 0; k < particleCount; k++)
                {
                    if(k==n) continue;
                    for (int i = 0; i < 2; i++)
                    {
                        tmp1 = roundParticles[n].ParticleDistance(roundParticles[k]);
                        if(tmp1 < particleViewRange)
                        {
                            roundParticles[n].forceFromCollisions[i] += (roundParticles[n].pos[i] - roundParticles[k].pos[i])*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
                            averageMuki += new Vector2(Mathf.Sin(roundParticles[k].theta),Mathf.Cos(roundParticles[k].theta));
                            closeParticleCount++;
                        }
                    }
                }
                if(closeParticleCount!=0)
                {
                    averageMuki /= closeParticleCount;
                    roundParticles[n].theta = Mathf.Atan2(averageMuki[0],averageMuki[1]) + (Random.Range(-noiseAngleDeg,noiseAngleDeg)*Mathf.PI)/180f;
                }
            }
            else
            {
                for (int k = 0; k < particleCount; k++)
                {
                    if(k==n) continue;
                    for (int i = 0; i < 2; i++)
                    {
                        tmp1 = roundParticles[n].ParticleDistance(roundParticles[k]);
                        if(tmp1 < particleViewRange)
                        {
                            roundParticles[n].forceFromCollisions[i] += (roundParticles[n].pos[i] - roundParticles[k].pos[i])*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
                        }
                    }
                }
            }
            

            float wallDistance;
            float angle;
            wallDistance = Mathf.Abs(roundParticles[n].pos[1]);
            if(wallDistance < particleViewRange)
            {
                if(wallDistance < roundParticles[n].radius) roundParticles[n].pos[1] = roundParticles[n].radius;
                angle = Vector2.SignedAngle(
                    new Vector2(0,-1),
                    new Vector2(-Mathf.Sin(roundParticles[n].theta),Mathf.Cos(roundParticles[n].theta))
                );
                if(Mathf.Abs(angle) < particleFov/2f)
                {
                    float turnAngle = particleFov/2f - Mathf.Abs(angle);
                    if(angle < 0) turnAngle *= -1;
                    roundParticles[n].theta += (turnAngle*Mathf.PI)/180f;
                }
            }
            wallDistance = Mathf.Abs(DIM_Y-1-roundParticles[n].pos[1]);
            if(wallDistance < particleViewRange)
            {
                if(wallDistance < roundParticles[n].radius) roundParticles[n].pos[1] = DIM_Y-1-roundParticles[n].radius;
                angle = Vector2.SignedAngle(
                    new Vector2(0,1),
                    new Vector2(-Mathf.Sin(roundParticles[n].theta),Mathf.Cos(roundParticles[n].theta))
                );
                if(Mathf.Abs(angle) < particleFov/2f)
                {
                    float turnAngle = particleFov/2f - Mathf.Abs(angle);
                    if(angle < 0) turnAngle *= -1;
                    roundParticles[n].theta +=(turnAngle*Mathf.PI)/180f;
                }
            }
            wallDistance = Mathf.Abs(roundParticles[n].pos[0]);
            if(wallDistance < particleViewRange)
            {
                if(wallDistance < roundParticles[n].radius) roundParticles[n].pos[0] =roundParticles[n].radius;
                angle = Vector2.SignedAngle(
                    new Vector2(-1,0),
                    new Vector2(-Mathf.Sin(roundParticles[n].theta),Mathf.Cos(roundParticles[n].theta))
                );
                if(Mathf.Abs(angle) < particleFov/2f)
                {
                    float turnAngle = particleFov/2f - Mathf.Abs(angle);
                    if(angle < 0) turnAngle *= -1;
                    roundParticles[n].theta += (turnAngle*Mathf.PI)/180f;
                }
            }
            wallDistance = Mathf.Abs(DIM_X-1-roundParticles[n].pos[0]);
            if(wallDistance < particleViewRange)
            {
                if(wallDistance < roundParticles[n].radius) roundParticles[n].pos[0] = DIM_X-1-roundParticles[n].radius;
                angle = Vector2.SignedAngle(
                    new Vector2(1,0),
                    new Vector2(-Mathf.Sin(roundParticles[n].theta),Mathf.Cos(roundParticles[n].theta))
                );
                if(Mathf.Abs(angle) < particleFov/2f)
                {
                    float turnAngle = particleFov/2f - Mathf.Abs(angle);
                    if(angle < 0) turnAngle *= -1;
                    roundParticles[n].theta +=(turnAngle*Mathf.PI)/180f;
                }
            }


            

            roundParticles[n].forceFromFluid[0] = 0f;
            roundParticles[n].forceFromFluid[1] = 0f;
            roundParticles[n].torque = 0f;
            for(int m = 0; m < roundParticles[n].perimeterPointCount ; m++) 
            {
                roundParticles[n].perimeterFluidVel[m,0] = 0f;
                roundParticles[n].perimeterFluidVel[m,1] = 0f;
                angle = Vector2.Angle(
                    new Vector2(roundParticles[n].perimeterPos[m,0] - roundParticles[n].pos[0],roundParticles[n].perimeterPos[m,1] - roundParticles[n].pos[1]),
                    new Vector2(-Mathf.Sin(roundParticles[n].theta),Mathf.Cos(roundParticles[n].theta))
                );
                // 固体表面の速度を計算
                for(int i = (int)roundParticles[n].perimeterPos[m,0] - 3; i < (int)roundParticles[n].perimeterPos[m,0] + 3; i++)
                {
                    for(int j = (int)roundParticles[n].perimeterPos[m,1] - 3; j < (int)roundParticles[n].perimeterPos[m,1] + 3; j++)
                    {
                        tmp1 = Mathf.Abs(roundParticles[n].perimeterPos[m,0] - (float)i);
                        tmp2 = Mathf.Abs(roundParticles[n].perimeterPos[m,1] - (float)j);
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
                            roundParticles[n].perimeterFluidVel[m,0] += u[i + j*DIM_X]*tmp3;
                            roundParticles[n].perimeterFluidVel[m,1] += v[i + j*DIM_X]*tmp3;
                        }
                    } 
                }
                roundParticles[n].forceOnPerimeter[m,0] = roundParticles[n].perimeterVel[m,0] - roundParticles[n].perimeterFluidVel[m,0];
                roundParticles[n].forceOnPerimeter[m,1] = roundParticles[n].perimeterVel[m,1] - roundParticles[n].perimeterFluidVel[m,1];

                if(angle>90f)
                {
                    Vector2 vect2Center = new Vector2(
                        roundParticles[n].pos[0] - roundParticles[n].perimeterPos[m,0],
                        roundParticles[n].pos[1] - roundParticles[n].perimeterPos[m,1]
                    );
                    vect2Center.Normalize();
                    roundParticles[n].forceOnPerimeter[m,0] += particlePropulsion * vect2Center[0] * Mathf.Cos((angle*Mathf.PI)/180f);
                    roundParticles[n].forceOnPerimeter[m,1] += particlePropulsion * vect2Center[1] * Mathf.Cos((angle*Mathf.PI)/180f);
                }

                // 固体が外部に与える力を計算
                for(int i = (int)roundParticles[n].perimeterPos[m,0] - 3; i < (int)roundParticles[n].perimeterPos[m,0] + 3; i++)
                {
                    for(int j = (int)roundParticles[n].perimeterPos[m,1] - 3; j < (int)roundParticles[n].perimeterPos[m,1] + 3; j++)
                    {
                        tmp1 = Mathf.Abs(roundParticles[n].perimeterPos[m,0] - (float)i);
                        tmp2 = Mathf.Abs(roundParticles[n].perimeterPos[m,1] - (float)j);
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
                            forceFromParticlesX[i + j*DIM_X] += roundParticles[n].forceOnPerimeter[m,0] * tmp3 * 2.0f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;
                            forceFromParticlesY[i + j*DIM_X] += roundParticles[n].forceOnPerimeter[m,1] * tmp3 * 2.0f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;
                        }
                    } 
                }
                roundParticles[n].forceFromFluid[0] += roundParticles[n].forceOnPerimeter[m,0];
                roundParticles[n].forceFromFluid[1] += roundParticles[n].forceOnPerimeter[m,1];
                roundParticles[n].torque += roundParticles[n].forceOnPerimeter[m,1] * (roundParticles[n].perimeterPos[m,0] - roundParticles[n].pos[0]) 
                                        - roundParticles[n].forceOnPerimeter[m,0] * (roundParticles[n].perimeterPos[m,1] - roundParticles[n].pos[1]);
            } 

            roundParticles[n].forceFromFluid[0] *= -2f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;  
            roundParticles[n].forceFromFluid[1] *= -2f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;  
            roundParticles[n].torque *= -2f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;  

            // roundParticles[n].UpdatePosVel(new float[2]{-Mathf.Sin(roundParticles[n].theta)*particlePropulsion,Mathf.Cos(roundParticles[n].theta)*particlePropulsion});
            roundParticles[n].UpdatePosVel();
            roundParticles[n].UpdateOmegaTheta();
            roundParticles[n].UpdatePerimeter();
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

    // public void Reset()
    // {
    //     DIM_X = int.Parse(inputFieldW.text);
    //     DIM_Y = int.Parse(inputFieldH.text);
    //     time = 0;
    //     plotTexture = new Texture2D(DIM_X,DIM_Y);
    //     plotTexture.filterMode = FilterMode.Point;
    //     plotPixels = plotTexture.GetPixels();
    //     plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),UnityEngine.Vector2.zero);

    //     if(DIM_X>DIM_Y)
    //     {
    //         plotImage.transform.GetComponent<RectTransform>().sizeDelta = new Vector2(
    //             plotImage.transform.GetComponent<RectTransform>().sizeDelta.x,
    //             (plotImage.transform.GetComponent<RectTransform>().sizeDelta.x*(float)DIM_Y)/(float)DIM_X
    //         );
    //     }
    //     if(DIM_X<DIM_Y)
    //     {
    //         plotImage.transform.GetComponent<RectTransform>().sizeDelta = new Vector2(
    //             (plotImage.transform.GetComponent<RectTransform>().sizeDelta.y*(float)DIM_X)/(float)DIM_Y,
    //             plotImage.transform.GetComponent<RectTransform>().sizeDelta.y
    //         );
    //     }
    //     vectorField.DestroyAll();
    //     vectorField.aspectRatio = (float)DIM_X/(float)DIM_Y;
    //     vectorField.Start();
    //     speed = new float[DIM_X*DIM_Y];
    //     u = new float[DIM_X*DIM_Y];
    //     v = new float[DIM_X*DIM_Y];
    //     fx = new float[DIM_X*DIM_Y];
    //     forceFromParticlesX = new float[DIM_X*DIM_Y];
    //     forceFromParticlesY = new float[DIM_X*DIM_Y];
    //     fy = new float[DIM_X*DIM_Y];
    //     rho = new float[DIM_X*DIM_Y];
    //     temp = new float[DIM_X*DIM_Y];
    //     f = new float[9*DIM_X*DIM_Y];
    //     f0 = new float[9*DIM_X*DIM_Y];
    //     ftmp = new float[9*DIM_X*DIM_Y];
    //     g = new float[5*DIM_X*DIM_Y];
    //     g0 = new float[5*DIM_X*DIM_Y];
    //     gtmp = new float[5*DIM_X*DIM_Y];
    //     plotSizeDelta = plotImage.transform.GetComponent<RectTransform>().sizeDelta;
    //     plotScale = plotImage.transform.localScale;
    //     if(!isFaller)
    //     {
    //         float[] spawnPos = new float[particleCount*2];
    //         GameObject[] objs = new GameObject[particleCount];
    //         particleRadiusSlider.value = particleRadius;
    //         float scale = (2f*particleRadius*plotImage.transform.GetComponent<RectTransform>().sizeDelta.x * plotImage.transform.localScale.x)/(DIM_X*particlePrefab.GetComponent<RectTransform>().sizeDelta.x);
    //         particleCountSlider.value = particleCount;
    //         for (int i = 0; i < particleCount; i++)
    //         {
    //             Destroy(roundParticles.objs[i]);
    //             spawnPos[i + 0*particleCount] = UnityEngine.Random.Range(particleRadius,DIM_X-particleRadius);
    //             spawnPos[i + 1*particleCount] = UnityEngine.Random.Range(particleRadius,DIM_Y-particleRadius);
    //             objs[i] = Instantiate(particlePrefab,particleParent);
    //             objs[i].GetComponent<RectTransform>().anchoredPosition = (-plotSizeDelta/2f  
    //                                                             + new Vector2((spawnPos[i + 0*particleCount]*plotSizeDelta.x)/DIM_X,(spawnPos[i + 1*particleCount]*plotSizeDelta.y)/DIM_Y))*plotScale.x;
    //             objs[i].transform.localScale = new Vector3(scale,scale,1);
    //         }
    //         roundParticles = new RoundParticlesForCompute(DIM_X,particleCount,particleDensity,particleRadius,spawnPos,objs,Random.Range(0f,2f*Mathf.PI));
    //     }
    //     else
    //     {
    //         UpdatePolygonPerimeter(true);
    //     }

    //     maxTemp = 0f;
    //     minTemp = Mathf.Infinity;
    //     for(int i = 0; i < DIM_X; i++)
    //     { 
    //         for(int j = 0; j < DIM_Y; j++)
    //         {
    //             u[i + j*DIM_X] = 0.0f; v[i + j*DIM_X] = 0.0f; 
    //             rho[i + j*DIM_X] = 1.0f;

    //             u2 = u[i + j*DIM_X]*u[i + j*DIM_X] + v[i + j*DIM_X]*v[i + j*DIM_X];   

    //             forceFromParticlesX[i + j*DIM_X] = 0f;
    //             forceFromParticlesY[i + j*DIM_X] = 0f;
    //             if(initModeDropDown.value == 0) temp[i + j*DIM_X] = 0f;
    //             else if(initModeDropDown.value == 1) temp[i + j*DIM_X] = (float)(DIM_X-1 - i)/(float)(DIM_X-1 - 1);
    //             else if(initModeDropDown.value == 2) temp[i + j*DIM_X] = (float)(DIM_Y-1 - j)/(float)(DIM_Y-1 - 1);
                
    //             for (int k = 0; k < 9; k++)
    //             {
    //                 tmp = cx[k]*u[i + j*DIM_X] + cy[k]*v[i + j*DIM_X];     
    //                 f0[k + (i + j*DIM_X)*9] = wf[k]*rho[i + j*DIM_X]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
    //                 f[k + (i + j*DIM_X)*9] = f0[k + (i + j*DIM_X)*9];
    //                 if(k<5)
    //                 {
    //                     g0[k + (i + j*DIM_X)*5] = wg[k]*temp[i + j*DIM_X]*(1.0f + 3.0f*tmp);
    //                     g[k + (i + j*DIM_X)*5] = g0[k + (i + j*DIM_X)*5];
    //                 }
    //             }
    //         }
    //     }
    // }

    public void ZeroGrav()
    {
        gxSlider.value = 0f;
        gySlider.value = 0f;
    }
}