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

public class RoundParticle
{
    public float density;
    public float radius;
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
    public RoundParticle(float _density,float _radius,float[] _initPos)
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
            // // perimeterVel[i,0] = vel[0] - omega*(perimeterPos[i,1] - pos[1]);
            // // perimeterVel[i,1] = vel[1] + omega*(perimeterPos[i,0] - pos[0]);
            perimeterVel[i,0] = 0f;
            perimeterVel[i,1] = 0f;
            forceOnPerimeter[i,0] = 0f;
            forceOnPerimeter[i,1] = 0f;
            perimeterFluidVel[i,0] = 0f;
            perimeterFluidVel[i,1] = 0f;
        } 
    }

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
}
