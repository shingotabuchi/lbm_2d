using UnityEngine;
using UnityEngine.UI;

public class VectorField : MonoBehaviour
{
    public RectTransform plotImageRect;
    public GameObject vectorImage;
    public Transform canvasTrans;
    public int resInt;
    int resIntX;
    public float scale = 1f;
    public GameObject[,] vectors;
    public float aspectRatio = 1f;
    float dx,dy,rangeX,rangeY;
    void Start()
    {
        resIntX = (int)(resInt*aspectRatio);
        vectors = new GameObject[resIntX+1,resInt+1];
        rangeX = plotImageRect.localScale.x * plotImageRect.sizeDelta.x;
        rangeY = plotImageRect.localScale.y * plotImageRect.sizeDelta.y;
        dx = rangeX/resIntX;
        dy = rangeY/resInt;
        Vector2 origin = plotImageRect.anchoredPosition - new Vector2(rangeX/2, rangeY/2);
        for (int i = 0; i <= resIntX; i++)
        {
            for (int j = 0; j <= resInt; j++)
            {
                vectors[i,j] = Instantiate(vectorImage,canvasTrans);
                vectors[i,j].GetComponent<RectTransform>().anchoredPosition = new Vector2(dx*i,dy*j) + origin;
                vectors[i,j].GetComponent<RectTransform>().sizeDelta = new Vector2(0,0);
            }
        }
    }

    public void UpdateVectors(float[,] u, float[,] v, int plotRes,float maxSpeed, bool on = true)
    {
        if(!on){
            canvasTrans.transform.gameObject.SetActive(false);
            return;
        }
        else canvasTrans.transform.gameObject.SetActive(true);
        for (int i = 0; i <= resIntX; i++)
        {
            for (int j = 0; j <= resInt; j++)
            {
                float U = u[(i*(int)(plotRes*aspectRatio-1))/resIntX,(j*(plotRes-1))/resInt];
                float V = v[(i*(int)(plotRes*aspectRatio-1))/resIntX,(j*(plotRes-1))/resInt];
                if(!(U==0f&&V==0f))
                {
                    float rad = Mathf.Atan(V/U);
                    if(U<0) rad += Mathf.PI;
                    vectors[i,j].transform.rotation = Quaternion.Euler(0,0,(rad*180f)/Mathf.PI);
                    float speed = Mathf.Sqrt(U*U + V*V);
                    if(speed>Mathf.Epsilon)
                    {
                        vectors[i,j].GetComponent<RectTransform>().sizeDelta = new Vector2(1200,1200) * (speed/maxSpeed) * scale;
                    }
                    else vectors[i,j].GetComponent<RectTransform>().sizeDelta = new Vector2(0,0);
                }
                else vectors[i,j].GetComponent<RectTransform>().sizeDelta = new Vector2(0,0);
            }
        }
    }
    public void UpdateVectors(float[] u, float[] v, int plotRes,float maxSpeed, bool on = true)
    {
        if(!on){
            canvasTrans.transform.gameObject.SetActive(false);
            return;
        }
        else canvasTrans.transform.gameObject.SetActive(true);
        for (int i = 0; i <= resIntX; i++)
        {
            for (int j = 0; j <= resInt; j++)
            {
                // float U = u[(i*(int)(plotRes*aspectRatio-1))/resIntX + ((j*(plotRes-1))/resInt)*resIntX];
                // float V = v[(i*(int)(plotRes*aspectRatio-1))/resIntX + ((j*(plotRes-1))/resInt)*resIntX];
                float U = u[(i*(int)(plotRes*aspectRatio-1))/resIntX + ((j*(plotRes-1))/resInt)*(int)(plotRes*aspectRatio)];
                float V = v[(i*(int)(plotRes*aspectRatio-1))/resIntX + ((j*(plotRes-1))/resInt)*(int)(plotRes*aspectRatio)];
                if(!(U==0f&&V==0f))
                {
                    float rad = Mathf.Atan(V/U);
                    if(U<0) rad += Mathf.PI;
                    vectors[i,j].transform.rotation = Quaternion.Euler(0,0,(rad*180f)/Mathf.PI);
                    float speed = Mathf.Sqrt(U*U + V*V);
                    if(speed>Mathf.Epsilon)
                    {
                        vectors[i,j].GetComponent<RectTransform>().sizeDelta = new Vector2(1200,1200) * (speed/maxSpeed) * scale;
                    }
                    else vectors[i,j].GetComponent<RectTransform>().sizeDelta = new Vector2(0,0);
                }
                else vectors[i,j].GetComponent<RectTransform>().sizeDelta = new Vector2(0,0);
            }
        }
    }
}