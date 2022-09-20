using UnityEngine;

public enum BoundaryType
{
    Constant,
    Adiabatic,
    Bounceback,
}

public enum BoundaryTypeTherm
{
    Dirichlet,
    Neumann,
}

public enum HeatMapMode
{
    Speed,
    Density,
    Temperature,
    ChemicalPotential,
    OrderParameter,
}

public enum RTType
{
    SRT,
    MRT,
}

public enum InitialTemperatureDistribution
{
    AllZero,
    XGradient,
    YGradient,
}

public class RectObstacle
{
    public int[] pos = new int[2];
    public int h,w;
    public int dim;
    public RectObstacle(int[] _pos,int _h, int _w,int _dim)
    {
        pos = _pos;
        h = _h;
        w = _w;
        dim = _dim;
    }

    public bool PointIsInRect(int x,int y)
    {
        int leftx = (pos[0]-w/2+dim)%dim;
        int rightx = (pos[0]+w/2+dim)%dim;

        if(leftx<=rightx){
            if(x<leftx) return false;
            if(x>rightx) return false;
        }
        else
        {
            if(x>rightx&&x<leftx) return false;
        }
        if(y>pos[1]+h/2) return false;
        if(y<pos[1]-h/2) return false;
        return true;
    }
}
public class ParticlesForCompute
{
    public int particleCount;
    public int DIM_X;
    public GameObject[] objs;
    public float[] density;
    public float[] omega;//角速度
    public float[] theta;//角
    public float[] prevOmega1;//前フレームの角速度
    public float[] prevOmega2;//前々フレームの角速度
    public float[] torque;
    public float[] pos = new float[2];
    public float[] vel = new float[2];
    public float[] prevVel1 = new float[2];//前フレームの速度
    public float[] prevVel2 = new float[2];//前々フレームの速度
    public float[] forceFromCollisions = new float[2];
    public float[] forceFromFluid = new float[2];
    public float[] perimeterPos;
    public float[] perimeterVel;
    public float[] perimeterFluidVel;
    public float[] forceOnPerimeter;
    public int[] perimeterPointCount;
    public float[] volume;
    public float[] mass;
    public float[] momentOfInertia;

    public void PlotParticlePerimeter(ref Color[] pixels, Color? color = null)
    {
        color ??= Color.white;
        for (int n = 0; n < particleCount; n++)
        {
            for(int i = 0; i < perimeterPointCount[n]; i++) 
            {
                if(
                    (int)perimeterPos[n + (i + 0*DIM_X)*particleCount] + (int)perimeterPos[n + (i + 1*DIM_X)*particleCount] * DIM_X < pixels.Length
                    &&  
                    (int)perimeterPos[n + (i + 0*DIM_X)*particleCount] + (int)perimeterPos[n + (i + 1*DIM_X)*particleCount] * DIM_X >= 0
                )
                pixels[(int)perimeterPos[n + (i + 0*DIM_X)*particleCount] + (int)perimeterPos[n + (i + 1*DIM_X)*particleCount] * DIM_X] = (Color)color;
            }
        }
    }
    public void PlotParticleTrace(ref Color[] pixels, Color? color = null)
    {
        color ??= Color.white;
        for (int n = 0; n < particleCount; n++)
        {
            for(int i = 0; i < perimeterPointCount[n] ; i++) 
            {
                if(
                    (int)perimeterPos[n + (i + 0*DIM_X)*particleCount] + (int)perimeterPos[n + (i + 1*DIM_X)*particleCount] * DIM_X < pixels.Length
                    &&  
                    (int)perimeterPos[n + (i + 0*DIM_X)*particleCount] + (int)perimeterPos[n + (i + 1*DIM_X)*particleCount] * DIM_X >= 0
                )
                pixels[(int)perimeterPos[n + (i + 0*DIM_X)*particleCount] + (int)perimeterPos[n + (i + 1*DIM_X)*particleCount] * DIM_X] = (Color)color;
            }
        }
    }
    public void UpdatePosVel(float[]? gravity = null)
    {
        gravity ??= new float[2]{0f,0f};
        for (int n = 0; n < particleCount; n++)
        {
            for (int i = 0; i < 2; i++)
            {
                vel[n + i*particleCount] = (1f + 1f/density[n]) * prevVel1[n + i*particleCount]
                                        - 1f/density[n] * prevVel2[n + i*particleCount]
                                        + (forceFromFluid[n + i*particleCount] + forceFromCollisions[n + i*particleCount])/mass[n]
                                        + (1f - 1f/density[n]) * (float)gravity[i];
                pos[n + i*particleCount] += (vel[n + i*particleCount] + prevVel1[n + i*particleCount])/2f;
                prevVel2[n + i*particleCount] = prevVel1[n + i*particleCount];
                prevVel1[n + i*particleCount] = vel[n + i*particleCount];
            }
        }
    }
    public void UpdateOmegaTheta()
    {
        for (int n = 0; n < particleCount; n++)
        {
            omega[n] = (1f + 1f/density[n]) * prevOmega1[n] 
                                    - 1f/density[n] * prevOmega2[n]
                                    + torque[n]/momentOfInertia[n];
            theta[n] += (omega[n] - prevOmega1[n])/2f;
            prevOmega2[n] = prevOmega1[n];
            prevOmega1[n] = omega[n];
        } 
    }
}
public class Particle
{
    public GameObject obj;
    public float density;
    public float omega;//角速度
    public float theta;//角
    public float prevOmega1;//前フレームの角速度
    public float prevOmega2;//前々フレームの角速度
    public float torque;
    public float[] pos = new float[2];
    public float[] vel = new float[2];
    public float[] prevVel1 = new float[2];//前フレームの速度
    public float[] prevVel2 = new float[2];//前々フレームの速度
    public float[] forceFromCollisions = new float[2];
    public float[] forceFromFluid = new float[2];
    public float[,] perimeterPos;
    public float[,] perimeterVel;
    public float[,] perimeterFluidVel;
    public float[,] forceOnPerimeter;
    public int perimeterPointCount;
    public float volume;
    public float mass;
    public float momentOfInertia;

    public void PlotParticlePerimeter(ref Color[] pixels, int dim_x, Color? color = null)
    {
        color ??= Color.white;
        for(int i = 0; i < perimeterPointCount ; i++) 
        {
            if(
                (int)perimeterPos[i,0] + (int)perimeterPos[i,1] * dim_x < pixels.Length
                &&  
                (int)perimeterPos[i,0] + (int)perimeterPos[i,1] * dim_x >= 0
            )
            pixels[(int)perimeterPos[i,0] + (int)perimeterPos[i,1] * dim_x] = (Color)color;
        }
    }
    public void PlotParticleTrace(ref Color[] pixels, int dim_x, Color? color = null)
    {
        color ??= Color.white;
        for(int i = 0; i < perimeterPointCount ; i++) 
        {
            if(
                (int)perimeterPos[i,0] + (int)perimeterPos[i,1] * dim_x < pixels.Length
                &&  
                (int)perimeterPos[i,0] + (int)perimeterPos[i,1] * dim_x >= 0
            )
            pixels[(int)perimeterPos[i,0] + (int)perimeterPos[i,1] * dim_x] = (Color)color;
        }
    }
        public void UpdatePosVel(float[]? gravity = null)
    {
        gravity ??= new float[2]{0f,0f};
        for (int i = 0; i < 2; i++)
        {
            vel[i] = (1f + 1f/density) * prevVel1[i]
                                    - 1f/density * prevVel2[i]
                                    + (forceFromFluid[i] + forceFromCollisions[i])/mass
                                    + (1f - 1f/density) * (float)gravity[i];
            pos[i] += (vel[i] + prevVel1[i])/2f;
            prevVel2[i] = prevVel1[i];
            prevVel1[i] = vel[i];
        }
    }
    public void UpdateOmegaTheta()
    {
        omega = (1f + 1f/density) * prevOmega1 
                                - 1f/density * prevOmega2
                                + torque/momentOfInertia;
        theta += (omega - prevOmega1)/2f;
        prevOmega2 = prevOmega1;
        prevOmega1 = omega;
    }

}

public class PolygonParticle : Particle
{
    public PolygonParticle(int _perimeterPointCount,float _density = 1f, float _area = 100f,float[]? _initPos = null)
    {
        perimeterPointCount = _perimeterPointCount;
        vel = new float[2]{0f,0f};
        prevVel1 = new float[2]{0f,0f};
        prevVel2 = new float[2]{0f,0f};
        omega = 0f;
        theta = 0f;
        prevOmega1 = 0f;
        prevOmega2 = 0f;
        perimeterPos = new float[perimeterPointCount,2];
        perimeterVel = new float[perimeterPointCount,2];
        perimeterFluidVel = new float[perimeterPointCount,2];
        forceOnPerimeter = new float[perimeterPointCount,2];
        forceFromCollisions = new float[2]{0f,0f};
        forceFromFluid = new float[2]{0f,0f};
        torque = 0f;
        volume = _area;
        density = _density;
        mass = volume * density;
        // momentOfInertia = (volume*volume/Mathf.PI)* _density/2f;
        momentOfInertia = 5100f;
        if(_initPos==null) pos = new float[2]{0f,0f};
        else pos = (float[])_initPos;
        for (int i = 0; i < _perimeterPointCount; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                perimeterVel[i,j] = 0f;
                perimeterFluidVel[i,j] = 0f;
                forceOnPerimeter[i,j] = 0f;
            }
        }
    }
    public PolygonParticle(int _perimeterPointCount, float[,] _perimeterPos, float _density, float _area)
    {
        perimeterPointCount = _perimeterPointCount;
        perimeterPos = _perimeterPos;
        volume = _area;
        density = _density;
        mass = volume * density;
        momentOfInertia = (volume*volume/Mathf.PI)* _density/2f;
    }
    public PolygonParticle()
    {
    }
}
public class RoundParticlesForCompute : ParticlesForCompute
{
    public float[] radius;
    public int maxPerimeterPointCount = 0;
    // public RoundParticlesForCompute(int _particleCount,float _density,float _radius,float[] _initPos, float _initTheta = 0f)
    public RoundParticlesForCompute(int _particleCount,float _density,float _radius,float[] _initPos,GameObject[] _objs, float _initTheta = 0f)
    {
        particleCount = _particleCount;
        density = new float[particleCount];
        radius = new float[particleCount];
        volume = new float[particleCount];
        mass = new float[particleCount];
        momentOfInertia = new float[particleCount];
        pos = new float[particleCount*2];
        vel = new float[particleCount*2];
        prevVel1 = new float[particleCount*2];
        prevVel2 = new float[particleCount*2];
        omega = new float[particleCount];
        theta = new float[particleCount];
        prevOmega1 = new float[particleCount];
        prevOmega2 = new float[particleCount];
        perimeterPointCount = new int[particleCount];
        forceFromCollisions = new float[particleCount*2];
        forceFromFluid = new float[particleCount*2];
        torque = new float[particleCount];
        pos = _initPos;
        objs = _objs;
        for (int n = 0; n < particleCount; n++)
        {
            density[n] = _density;
            radius[n] = _radius;
            volume[n] = Mathf.PI*radius[n]*radius[n];
            mass[n] = volume[n] * density[n];
            momentOfInertia[n] = (Mathf.PI * radius[n]*radius[n]*radius[n]*radius[n] * density[n])/2f;
            // pos[n + 0*particleCount] = _initPos[n + 0*particleCount];
            // pos[n + 1*particleCount] = _initPos[n + 1*particleCount];
            vel[n + 0*particleCount] = 0f;
            vel[n + 1*particleCount] = 0f;
            prevVel1[n + 0*particleCount] = 0f;
            prevVel1[n + 1*particleCount] = 0f;
            prevVel2[n + 0*particleCount] = 0f;
            prevVel2[n + 1*particleCount] = 0f;
            omega[n] = 0f;
            theta[n] = _initTheta;
            prevOmega1[n] = 0f;
            prevOmega2[n] = 0f;
            perimeterPointCount[n] = (int)(2.0f * Mathf.PI * radius[n] * 2.0f);
            maxPerimeterPointCount = (int)Mathf.Abs(perimeterPointCount[n]);
            forceFromCollisions[n + 0*particleCount] = 0f;
            forceFromCollisions[n + 1*particleCount] = 0f;
            forceFromFluid[n + 0*particleCount] = 0f;
            forceFromFluid[n + 1*particleCount] = 0f;
            torque[n] = 0f;
        }
        perimeterPos = new float[maxPerimeterPointCount*2*particleCount];
        perimeterVel = new float[maxPerimeterPointCount*2*particleCount];
        perimeterFluidVel = new float[maxPerimeterPointCount*2*particleCount];
        forceOnPerimeter = new float[maxPerimeterPointCount*2*particleCount];
        for (int n = 0; n < particleCount; n++)
        {
            for (int i = 0; i < perimeterPointCount[n]; i++)
            {
                float angle = (3f*Mathf.PI)/2f + 2.0f*Mathf.PI*(float)i/(float)perimeterPointCount[n] + theta[n];

                perimeterPos[n + (i + 0*maxPerimeterPointCount)*particleCount] = pos[n + 0*particleCount] + radius[n] * Mathf.Cos(angle);
                perimeterPos[n + (i + 1*maxPerimeterPointCount)*particleCount] = pos[n + 1*particleCount] + radius[n] * Mathf.Sin(angle);
                perimeterVel[n + (i + 0*maxPerimeterPointCount)*particleCount] = vel[n + 0*particleCount] - omega[n]*(perimeterPos[n + (i + 1*maxPerimeterPointCount)*particleCount] - pos[n + 1*particleCount]);
                perimeterVel[n + (i + 1*maxPerimeterPointCount)*particleCount] = vel[n + 1*particleCount] + omega[n]*(perimeterPos[n + (i + 0*maxPerimeterPointCount)*particleCount] - pos[n + 0*particleCount]);
                forceOnPerimeter[n + (i + 0*maxPerimeterPointCount)*particleCount] = 0f;
                forceOnPerimeter[n + (i + 1*maxPerimeterPointCount)*particleCount] = 0f;
                perimeterFluidVel[n + (i + 0*maxPerimeterPointCount)*particleCount] = 0f;
                perimeterFluidVel[n + (i + 1*maxPerimeterPointCount)*particleCount] = 0f;
            }
        }
    }

    
    public void PlotParticleFill(ref Color[] pixels, Color? color = null)
    {
        color ??= Color.white;
        for (int n = 0; n < particleCount; n++)
        {
            for (int i = -((int)radius[n]); i <= (int)radius[n]; i++)
            {
                for (int j = -((int)radius[n]); j <= (int)radius[n]; j++)
                {
                    if(i*i + j*j <= radius[n]*radius[n])
                    {
                        if(
                            i + (int)pos[n + 0*particleCount] + (j + (int)pos[n + 1*particleCount])* DIM_X < pixels.Length
                            &&  
                            i + (int)pos[n + 0*particleCount] + (j + (int)pos[n + 1*particleCount])* DIM_X >= 0
                        )
                        pixels[i + (int)pos[n + 0*particleCount] + (j + (int)pos[n + 1*particleCount])* DIM_X] = (Color)color;
                    }
                }
            }
        }   
    }

    public void UpdatePerimeter()
    {
        for (int n = 0; n < particleCount; n++)
        {
            for(int i = 0; i < perimeterPointCount[n]; i++) 
            {
                float angle = (3f*Mathf.PI)/2f + 2.0f*Mathf.PI*(float)i/(float)perimeterPointCount[n] + theta[n];
                perimeterPos[n + (i + 0*maxPerimeterPointCount)*particleCount] = pos[n + 0*particleCount] + radius[n] * Mathf.Cos(angle);
                perimeterPos[n + (i + 1*maxPerimeterPointCount)*particleCount] = pos[n + 1*particleCount] + radius[n] * Mathf.Sin(angle);
                perimeterVel[n + (i + 0*maxPerimeterPointCount)*particleCount] = vel[n + 0*particleCount] - omega[n]*(perimeterPos[n + (i + 1*maxPerimeterPointCount)*particleCount] - pos[n + 1*particleCount]);
                perimeterVel[n + (i + 1*maxPerimeterPointCount)*particleCount] = vel[n + 1*particleCount] + omega[n]*(perimeterPos[n + (i + 0*maxPerimeterPointCount)*particleCount] - pos[n + 0*particleCount]);
            } 
        }
    }

    public float ParticleDistance(int a, int b)
    {
        return Mathf.Sqrt( (pos[a + 0*particleCount]-pos[b + 0*particleCount])*(pos[a + 0*particleCount]-pos[b + 0*particleCount]) + (pos[a + 1*particleCount]-pos[b + 1*particleCount])*(pos[a + 1*particleCount]-pos[b + 1*particleCount]) );
    }

    public bool PointIsInParticle(int[] point)
    {
        for (int n = 0; n < particleCount; n++)
        {
            if(
                (point[n + 0*particleCount] - pos[n + 0*particleCount])*(point[n + 0*particleCount] - pos[n + 0*particleCount]) + (point[n + 1*particleCount] - pos[n + 1*particleCount])*(point[n + 1*particleCount] - pos[n + 1*particleCount])
                < radius[n]*radius[n]
            )return true;
        }
        return false;
    }
}
public class RoundParticle : Particle
{
    public float radius;
    public int bottomPointIndex = 0;
    public RoundParticle(float _density,float _radius,float[] _initPos, float _initTheta = 0f)
    {
        density = _density;
        radius = _radius;
        volume = Mathf.PI*radius*radius;
        mass = volume * density;
        momentOfInertia = (Mathf.PI * radius*radius*radius*radius * density)/2f;
        pos = _initPos;
        vel = new float[2]{0f,0f};
        prevVel1 = new float[2]{0f,0f};
        prevVel2 = new float[2]{0f,0f};
        omega = 0f;
        theta = _initTheta;
        prevOmega1 = 0f;
        prevOmega2 = 0f;
        perimeterPointCount = (int)(2.0f * Mathf.PI * radius * 2.0f);
        perimeterPos = new float[perimeterPointCount,2];
        perimeterVel = new float[perimeterPointCount,2];
        perimeterFluidVel = new float[perimeterPointCount,2];
        forceOnPerimeter = new float[perimeterPointCount,2];
        forceFromCollisions = new float[2]{0f,0f};
        forceFromFluid = new float[2]{0f,0f};
        torque = 0f;

        for(int i = 0; i < perimeterPointCount; i++) 
        {
            float angle = (3f*Mathf.PI)/2f + 2.0f*Mathf.PI*(float)i/(float)perimeterPointCount + theta;

            perimeterPos[i,0] = pos[0] + radius * Mathf.Cos(angle);
            perimeterPos[i,1] = pos[1] + radius * Mathf.Sin(angle);
            perimeterVel[i,0] = vel[0] - omega*(perimeterPos[i,1] - pos[1]);
            perimeterVel[i,1] = vel[1] + omega*(perimeterPos[i,0] - pos[0]);
            forceOnPerimeter[i,0] = 0f;
            forceOnPerimeter[i,1] = 0f;
            perimeterFluidVel[i,0] = 0f;
            perimeterFluidVel[i,1] = 0f;
        } 
    }

    
    public void PlotParticleFill(ref Color[] pixels, int dim_x, Color? color = null)
    {
        color ??= Color.white;
        for (int i = -((int)radius); i <= (int)radius; i++)
        {
            for (int j = -((int)radius); j <= (int)radius; j++)
            {
                if(i*i + j*j <= radius*radius)
                {
                    if(
                        i + (int)pos[0] + (j + (int)pos[1])* dim_x < pixels.Length
                        &&  
                        i + (int)pos[0] + (j + (int)pos[1])* dim_x >= 0
                    )
                    pixels[i + (int)pos[0] + (j + (int)pos[1])* dim_x] = (Color)color;
                }
            }
        }
    }

    public void UpdatePosVelPeriodicX(int DIM_X)
    {
        for (int i = 0; i < 2; i++)
        {
            vel[i] = (1f + 1f/density) * prevVel1[i]
                                    - 1f/density * prevVel2[i]
                                    + (forceFromFluid[i] + forceFromCollisions[i])/mass;
                                    // + (1f - 1f/density) * gravity[i];
            pos[i] += (vel[i] + prevVel1[i])/2f;
            if(i==0)
            {
                if(pos[i]>=DIM_X) pos[i]-=DIM_X;
                else if(pos[i]<0) pos[i]+=DIM_X;
            }
            prevVel2[i] = prevVel1[i];
            prevVel1[i] = vel[i];
        }
    }

    public void UpdatePerimeterPeriodicX(int DIM_X)
    {
        for(int i = 0; i < perimeterPointCount; i++) 
        {
            perimeterPos[i,0] = pos[0] + radius * Mathf.Cos(2.0f*Mathf.PI*(float)i/(float)perimeterPointCount);
            if(perimeterPos[i,0]>=DIM_X) perimeterPos[i,0]-=DIM_X;
            else if(perimeterPos[i,0]<0) perimeterPos[i,0]+=DIM_X;
            perimeterPos[i,1] = pos[1] + radius * Mathf.Sin(2.0f*Mathf.PI*(float)i/(float)perimeterPointCount);
            perimeterVel[i,0] = vel[0] - omega*(perimeterPos[i,1] - pos[1]);
            perimeterVel[i,1] = vel[1] + omega*(perimeterPos[i,0] - pos[0]);
        } 
    }

    public void UpdatePerimeter()
    {
        for(int i = 0; i < perimeterPointCount; i++) 
        {
            float angle = (3f*Mathf.PI)/2f + 2.0f*Mathf.PI*(float)i/(float)perimeterPointCount + theta;;
            perimeterPos[i,0] = pos[0] + radius * Mathf.Cos(angle);
            perimeterPos[i,1] = pos[1] + radius * Mathf.Sin(angle);
            perimeterVel[i,0] = vel[0] - omega*(perimeterPos[i,1] - pos[1]);
            perimeterVel[i,1] = vel[1] + omega*(perimeterPos[i,0] - pos[0]);
        } 
    }

    public float ParticleDistance(RoundParticle particle)
    {
        return Mathf.Sqrt( (pos[0]-particle.pos[0])*(pos[0]-particle.pos[0]) + (pos[1]-particle.pos[1])*(pos[1]-particle.pos[1]) );
    }

    public bool PointIsInParticle(int[] point)
    {
        if(
            (point[0] - pos[0])*(point[0] - pos[0]) + (point[1] - pos[1])*(point[1] - pos[1])
            < radius*radius
        )return true;
        return false;
    }
}