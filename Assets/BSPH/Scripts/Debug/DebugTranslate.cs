using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class DebugTranslate : MonoBehaviour
{
    public enum TranslationAxis {
        X,Y,Z
    }
    public TranslationAxis axis = TranslationAxis.Y;
    public float translationRate = 2f;
    public float minimumOffset = 1f, maximumOffset = 1f;
    private float min, max;

    // Start is called before the first frame update
    void Start()
    {
        switch(axis) {
            case TranslationAxis.X:
                min = transform.position.x - minimumOffset;
                max = transform.position.x + maximumOffset;
                break;
            case TranslationAxis.Z:
                min = transform.position.z - minimumOffset;
                max = transform.position.z + maximumOffset;
                break;
            default:
                min = transform.position.y - minimumOffset;
                max = transform.position.y + maximumOffset;
                break;
        }
    }

    // Update is called once per frame
    void Update()
    {
        switch(axis) {
            case TranslationAxis.X:
                TranslateX();
                break;
            case TranslationAxis.Z:
                TranslateZ();
                break;
            default:
                TranslateY();
                break;
        }
    }

    void TranslateX() {
        transform.position = new Vector3(
            Mathf.PingPong(Time.time*translationRate,max-min)+min,
            transform.position.y,
            transform.position.z
        );
    }
    void TranslateY() {
        transform.position = new Vector3(
            transform.position.x,
            Mathf.PingPong(Time.time*translationRate,max-min)+min,  
            transform.position.z
        );
    }
    void TranslateZ() {
        transform.position = new Vector3(
            transform.position.x,
            transform.position.y,
            Mathf.PingPong(Time.time*translationRate,max-min)+min
        );
    }
}
