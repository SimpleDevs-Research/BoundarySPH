using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FPS : MonoBehaviour
{
    private float refreshTime = 0.5f;
    private int frameCounter = 0;
    private float timeCounter = 0f;
    [SerializeField] private float lastFrameRate = 0f;

    // Update is called once per frame
    void Update() {
        if( timeCounter < refreshTime ) {
            timeCounter += Time.deltaTime;
            frameCounter++;
        }
        else {
            //This code will break if you set your m_refreshTime to 0, which makes no sense.
            lastFrameRate = (float)frameCounter/timeCounter;
            frameCounter = 0;
            timeCounter = 0.0f;
        }
    }
}
