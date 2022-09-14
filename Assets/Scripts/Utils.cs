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

public class Particle
{
    public GameObject obj;
    public float density;
    public float omega;//角速度
    public float theta;//角
    public float prevOmega1;//前フレームの角速度
    public float prevOmega2;//前々フレームの角速度
    public float torque;
    public float[] pos;
    public float[] vel;
    public float[] prevVel1;//前フレームの速度
    public float[] prevVel2;//前々フレームの速度
    public float[] forceFromCollisions;
    public float[] forceFromFluid;
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
    public PolygonParticle(int _perimeterPointCount, float[,] _perimeterPos)
    {
        perimeterPointCount = _perimeterPointCount;
        perimeterPos = _perimeterPos;
    }
}

public class RoundParticle : Particle
{
    public float radius;
    public RoundParticle(float _density,float _radius,float[] _initPos, float _initOmega = 0f)
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
        omega = _initOmega;
        theta = 0f;
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
            perimeterPos[i,0] = pos[0] + radius * Mathf.Cos(2.0f*Mathf.PI*(float)i/(float)perimeterPointCount);
            perimeterPos[i,1] = pos[1] + radius * Mathf.Sin(2.0f*Mathf.PI*(float)i/(float)perimeterPointCount);
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
            perimeterPos[i,0] = pos[0] + radius * Mathf.Cos(2.0f*Mathf.PI*(float)i/(float)perimeterPointCount);
            perimeterPos[i,1] = pos[1] + radius * Mathf.Sin(2.0f*Mathf.PI*(float)i/(float)perimeterPointCount);
            perimeterVel[i,0] = vel[0] - omega*(perimeterPos[i,1] - pos[1]);
            perimeterVel[i,1] = vel[1] + omega*(perimeterPos[i,0] - pos[0]);
        } 
    }

    public float ParticleDistance(RoundParticle particle)
    {
        return Mathf.Sqrt( (pos[0]-particle.pos[0])*(pos[0]-particle.pos[0]) + (pos[1]-particle.pos[1])*(pos[1]-particle.pos[1]) );
    }
}