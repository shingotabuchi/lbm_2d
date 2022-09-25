using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;

public class Menu : MonoBehaviour
{
    public void LoadScene(int index)
    {
        switch(index)
        {
            case 0:
            SceneManager.LoadScene("TC");
            break;
            
            case 1:
            SceneManager.LoadScene("NewFaller");
            break;

            case 2:
            SceneManager.LoadScene("VicsekLike");
            break;
        }
    }
}
